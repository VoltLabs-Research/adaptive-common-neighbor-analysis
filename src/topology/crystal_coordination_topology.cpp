#include <volt/topology/crystal_coordination_topology.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/cna_classifier.h>
#include <volt/analysis/crystal_symmetry_utils.h>
#include <volt/math/quaternion.h>

#include <limits>

namespace Volt{

bool orthonormalizeOrientation(const Matrix3& input, Matrix3& output){
    Vector3 c0 = input.column(0);
    Vector3 c1 = input.column(1);
    Vector3 c2 = input.column(2);

    const double l0 = c0.length();
    if(l0 <= EPSILON){
        return false;
    }
    c0 /= l0;

    c1 -= c0 * c0.dot(c1);
    const double l1 = c1.length();
    if(l1 <= EPSILON){
        return false;
    }
    c1 /= l1;

    c2 -= c0 * c0.dot(c2);
    c2 -= c1 * c1.dot(c2);
    const double l2 = c2.length();
    if(l2 <= EPSILON){
        return false;
    }
    c2 /= l2;

    output = Matrix3(c0, c1, c2);
    if(output.determinant() < 0.0){
        output.column(2) = -output.column(2);
    }
    return true;
}

bool latticeSupportsOrientationInversion(
    const LatticeStructure& lattice,
    int coordinationNumber,
    std::array<int, MAX_NEIGHBORS>& outInversePermutation
){
    outInversePermutation.fill(-1);
    for(int slot = 0; slot < coordinationNumber; ++slot){
        const Vector3& vector = lattice.latticeVectors[static_cast<std::size_t>(slot)];
        for(int candidateSlot = 0; candidateSlot < coordinationNumber; ++candidateSlot){
            if((lattice.latticeVectors[static_cast<std::size_t>(candidateSlot)] + vector).isZero(EPSILON)){
                outInversePermutation[static_cast<std::size_t>(slot)] = candidateSlot;
                break;
            }
        }
        if(outInversePermutation[static_cast<std::size_t>(slot)] < 0){
            return false;
        }
    }
    return true;
}

bool computeMatchedOrientation(
    const LatticeStructure& lattice,
    const Vector3* neighborVectors,
    const int* neighborMapping,
    int coordinationNumber,
    Matrix3& outRawOrientation,
    Matrix3& outOrientation
){
    Matrix3 orientationV = Matrix3::Zero();
    Matrix3 orientationW = Matrix3::Zero();
    int vectorCount = 0;

    for(int canonicalSlot = 0; canonicalSlot < coordinationNumber; ++canonicalSlot){
        const int localSlot = neighborMapping[canonicalSlot];
        if(localSlot < 0 || localSlot >= coordinationNumber){
            return false;
        }

        const Vector3& idealVector = lattice.latticeVectors[static_cast<std::size_t>(canonicalSlot)];
        const Vector3& spatialVector = neighborVectors[static_cast<std::size_t>(localSlot)];
        for(int row = 0; row < 3; ++row){
            for(int column = 0; column < 3; ++column){
                orientationV(row, column) += idealVector[column] * idealVector[row];
                orientationW(row, column) += idealVector[column] * spatialVector[row];
            }
        }
        ++vectorCount;
    }

    if(vectorCount < 3){
        return false;
    }

    Matrix3 orientationVInverse;
    if(!orientationV.inverse(orientationVInverse)){
        return false;
    }

    outRawOrientation = Matrix3(orientationW * orientationVInverse);
    return orthonormalizeOrientation(outRawOrientation, outOrientation);
}

bool canonicalizeNeighborMapping(
    const LatticeStructure& lattice,
    const Vector3* neighborVectors,
    int coordinationNumber,
    int* neighborMapping
){
    if(coordinationNumber <= 0 ||
       static_cast<int>(lattice.latticeVectors.size()) < coordinationNumber ||
       lattice.permutations.empty()){
        return true;
    }

    auto scoreOrientation = [](const Matrix3& orientation) {
        Quaternion q(orientation);
        q.normalize();
        if(q.w() < 0.0){
            q = -q;
        }
        return std::array<double, 4>{q.w(), q.x(), q.y(), q.z()};
    };

    auto isBetterScore = [](const std::array<double, 4>& lhs, const std::array<double, 4>& rhs) {
        constexpr double tolerance = 1e-12;
        for(std::size_t index = 0; index < lhs.size(); ++index){
            if(lhs[index] > rhs[index] + tolerance){
                return true;
            }
            if(lhs[index] + tolerance < rhs[index]){
                return false;
            }
        }
        return false;
    };

    std::array<int, MAX_NEIGHBORS> baseMapping{};
    std::copy_n(neighborMapping, coordinationNumber, baseMapping.begin());

    std::array<int, MAX_NEIGHBORS> inversePermutation{};
    const bool allowInversionEquivalent = latticeSupportsOrientationInversion(
        lattice,
        coordinationNumber,
        inversePermutation
    );

    std::array<int, MAX_NEIGHBORS> bestMapping{};
    std::array<double, 4> bestScore{};
    bool bestValid = false;
    bool bestProper = false;

    auto tryCandidate = [&](const std::array<int, MAX_NEIGHBORS>& canonicalPermutation) {
        std::array<int, MAX_NEIGHBORS> candidateMapping{};
        for(int canonicalSlot = 0; canonicalSlot < coordinationNumber; ++canonicalSlot){
            const int mappedCanonicalSlot = canonicalPermutation[static_cast<std::size_t>(canonicalSlot)];
            if(mappedCanonicalSlot < 0 || mappedCanonicalSlot >= coordinationNumber){
                return;
            }
            candidateMapping[static_cast<std::size_t>(canonicalSlot)] =
                baseMapping[static_cast<std::size_t>(mappedCanonicalSlot)];
        }

        Matrix3 rawOrientation;
        Matrix3 orientation;
        if(!computeMatchedOrientation(
            lattice,
            neighborVectors,
            candidateMapping.data(),
            coordinationNumber,
            rawOrientation,
            orientation
        )){
            return;
        }

        const bool properOrientation = rawOrientation.determinant() >= 0.0;
        const std::array<double, 4> score = scoreOrientation(orientation);
        if(!bestValid ||
           (properOrientation && !bestProper) ||
           (properOrientation == bestProper && isBetterScore(score, bestScore))){
            bestMapping = candidateMapping;
            bestScore = score;
            bestValid = true;
            bestProper = properOrientation;
        }
    };

    for(const auto& symmetry : lattice.permutations){
        tryCandidate(symmetry.permutation);
    }

    if(allowInversionEquivalent){
        for(const auto& symmetry : lattice.permutations){
            std::array<int, MAX_NEIGHBORS> invertedPermutation{};
            invertedPermutation.fill(-1);
            for(int canonicalSlot = 0; canonicalSlot < coordinationNumber; ++canonicalSlot){
                const int rotatedSlot = symmetry.permutation[static_cast<std::size_t>(canonicalSlot)];
                if(rotatedSlot < 0 || rotatedSlot >= coordinationNumber){
                    continue;
                }
                invertedPermutation[static_cast<std::size_t>(canonicalSlot)] =
                    inversePermutation[static_cast<std::size_t>(rotatedSlot)];
            }
            tryCandidate(invertedPermutation);
        }
    }

    if(!bestValid){
        return false;
    }

    std::copy_n(bestMapping.begin(), coordinationNumber, neighborMapping);
    return true;
}

CoordinationStructure CoordinationStructures::_coordinationStructures[NUM_COORD_TYPES];
LatticeStructure CoordinationStructures::_latticeStructures[NUM_LATTICE_TYPES];

CoordinationStructures::CoordinationStructures(
    ParticleProperty* structureTypes,
    LatticeStructureType inputCrystalType,
    bool identifyPlanarDefects,
    const SimulationCell& simCell
    ) : _structureTypes(structureTypes)
        , _inputCrystalType(inputCrystalType)
        , _identifyPlanarDefects(identifyPlanarDefects)
        , _simCell(simCell){}

double CoordinationStructures::determineLocalStructure(
    const NearestNeighborFinder& neighList,
    int particleIndex,
    int* outNeighborCount,
    int* outOrderedNeighborIndices
) const{
    assert(_structureTypes->getInt(particleIndex) == COORD_OTHER);
    CnaLocalStructureUtils::LocalStructureMatch match;
    if(!CnaLocalStructureUtils::determineLocalStructure(
        _inputCrystalType,
        _identifyPlanarDefects,
        neighList,
        particleIndex,
        _coordinationStructures,
        match
    )){
        return 0.0;
    }

    const LatticeStructure& lattice = getLatticeStruct(static_cast<int>(match.structureType));
    if(!canonicalizeNeighborMapping(
        lattice,
        match.neighborVectors.data(),
        match.coordinationNumber,
        match.neighborMapping.data()
    )){
        return 0.0;
    }

    _structureTypes->setInt(particleIndex, static_cast<int>(match.structureType));

    if(outNeighborCount){
        *outNeighborCount = match.coordinationNumber;
    }

    if(outOrderedNeighborIndices){
        for(int i = 0; i < match.coordinationNumber; ++i){
            outOrderedNeighborIndices[i] = match.neighborIndices[match.neighborMapping[i]];
        }
    }

    return match.localCutoff;
}

void CoordinationStructures::postProcessDiamondNeighbors(
    AnalysisContext& context,
    const NearestNeighborFinder& neighList
) const{
    if(_inputCrystalType != LATTICE_CUBIC_DIAMOND && _inputCrystalType != LATTICE_HEX_DIAMOND){
        return;
    }

    const size_t N = _structureTypes->size();
    std::vector<int> newStructureTypes(N);

    for(size_t i = 0; i < N; ++i){
        newStructureTypes[i] = _structureTypes->getInt(i);
    }

    NearestNeighborFinder firstNeighborFinder(4);
    firstNeighborFinder.prepare(context.positions, context.simCell);

    for(size_t i = 0; i < N; ++i){
        int currentType = _structureTypes->getInt(i);
        if(currentType != static_cast<int>(StructureType::CUBIC_DIAMOND) &&
           currentType != static_cast<int>(StructureType::HEX_DIAMOND)){
            continue;
        }

        NearestNeighborFinder::Query<4> query(firstNeighborFinder);
        query.findNeighbors(i, false);

        StructureType firstNeighType = currentType == static_cast<int>(StructureType::CUBIC_DIAMOND)
            ? StructureType::CUBIC_DIAMOND_FIRST_NEIGH
            : StructureType::HEX_DIAMOND_FIRST_NEIGH;

        for(const auto& neighbor : query.results()){
            if(_structureTypes->getInt(neighbor.index) == static_cast<int>(StructureType::OTHER)){
                newStructureTypes[neighbor.index] = static_cast<int>(firstNeighType);
            }
        }
    }

    for(size_t i = 0; i < N; ++i){
        _structureTypes->setInt(i, newStructureTypes[i]);
    }

    for(size_t i = 0; i < N; ++i){
        int currentType = _structureTypes->getInt(i);
        if(currentType != static_cast<int>(StructureType::CUBIC_DIAMOND_FIRST_NEIGH) &&
           currentType != static_cast<int>(StructureType::HEX_DIAMOND_FIRST_NEIGH)){
            continue;
        }

        NearestNeighborFinder::Query<4> query(firstNeighborFinder);
        query.findNeighbors(i, false);

        StructureType secondNeighType = currentType == static_cast<int>(StructureType::CUBIC_DIAMOND_FIRST_NEIGH)
            ? StructureType::CUBIC_DIAMOND_SECOND_NEIGH
            : StructureType::HEX_DIAMOND_SECOND_NEIGH;

        for(const auto& neighbor : query.results()){
            if(_structureTypes->getInt(neighbor.index) == static_cast<int>(StructureType::OTHER)){
                newStructureTypes[neighbor.index] = static_cast<int>(secondNeighType);
            }
        }
    }

    for(size_t i = 0; i < N; ++i){
        _structureTypes->setInt(i, newStructureTypes[i]);
    }
}

void CoordinationStructures::initializeFCC(){
    initializeCoordinationStructure(COORD_FCC, FCC_VECTORS, 12, [&](const Vector3& v1, const Vector3& v2){
        return (v1 - v2).length() < (std::sqrt(0.5) + 1.0) * 0.5;
    }, [](int){ return 0; });

    initializeLatticeStructure(LATTICE_FCC, FCC_VECTORS, 12, &_coordinationStructures[COORD_FCC]);
    _latticeStructures[LATTICE_FCC].primitiveCell.column(0) = FCC_PRIMITIVE_CELL[0];
    _latticeStructures[LATTICE_FCC].primitiveCell.column(1) = FCC_PRIMITIVE_CELL[1];
    _latticeStructures[LATTICE_FCC].primitiveCell.column(2) = FCC_PRIMITIVE_CELL[2];
}

void CoordinationStructures::initializeHCP(){
    initializeCoordinationStructure(COORD_HCP, HCP_VECTORS, 12, [&](const Vector3& v1, const Vector3& v2){
        return (v1 - v2).length() < (std::sqrt(0.5) + 1.0) * 0.5;
    }, [&](int ni){ return (HCP_VECTORS[ni].z() == 0) ? 1 : 0; });

    initializeLatticeStructure(LATTICE_HCP, HCP_VECTORS, 18, &_coordinationStructures[COORD_HCP]);
    _latticeStructures[LATTICE_HCP].primitiveCell.column(0) = HCP_PRIMITIVE_CELL[0];
    _latticeStructures[LATTICE_HCP].primitiveCell.column(1) = HCP_PRIMITIVE_CELL[1];
    _latticeStructures[LATTICE_HCP].primitiveCell.column(2) = HCP_PRIMITIVE_CELL[2];
}

void CoordinationStructures::initializeBCC(){
    initializeCoordinationStructure(COORD_BCC, BCC_VECTORS, 14, [&](const Vector3& v1, const Vector3& v2){
        return (v1 - v2).length() < (1.0 + std::sqrt(2.0)) * 0.5;
    }, [](int ni){ return (ni < 8) ? 0 : 1; });

    initializeLatticeStructure(LATTICE_BCC, BCC_VECTORS, 14, &_coordinationStructures[COORD_BCC]);
    _latticeStructures[LATTICE_BCC].primitiveCell.column(0) = BCC_PRIMITIVE_CELL[0];
    _latticeStructures[LATTICE_BCC].primitiveCell.column(1) = BCC_PRIMITIVE_CELL[1];
    _latticeStructures[LATTICE_BCC].primitiveCell.column(2) = BCC_PRIMITIVE_CELL[2];
}

void CoordinationStructures::initializeCubicDiamond(){
    initializeDiamondStructure(COORD_CUBIC_DIAMOND, LATTICE_CUBIC_DIAMOND, DIAMOND_CUBIC_VECTORS, 16, 20);

    _latticeStructures[LATTICE_CUBIC_DIAMOND].primitiveCell.column(0) = CUBIC_DIAMOND_PRIMITIVE_CELL[0];
    _latticeStructures[LATTICE_CUBIC_DIAMOND].primitiveCell.column(1) = CUBIC_DIAMOND_PRIMITIVE_CELL[1];
    _latticeStructures[LATTICE_CUBIC_DIAMOND].primitiveCell.column(2) = CUBIC_DIAMOND_PRIMITIVE_CELL[2];
}

void CoordinationStructures::initializeHexagonalDiamond(){
    initializeDiamondStructure(COORD_HEX_DIAMOND, LATTICE_HEX_DIAMOND, DIAMOND_HEX_VECTORS, 16, 32);

    _latticeStructures[LATTICE_HEX_DIAMOND].primitiveCell.column(0) = HEXAGONAL_DIAMOND_PRIMITIVE_CELL[0];
    _latticeStructures[LATTICE_HEX_DIAMOND].primitiveCell.column(1) = HEXAGONAL_DIAMOND_PRIMITIVE_CELL[1];
    _latticeStructures[LATTICE_HEX_DIAMOND].primitiveCell.column(2) = HEXAGONAL_DIAMOND_PRIMITIVE_CELL[2];
}

void CoordinationStructures::initializeSC(){
    initializeCoordinationStructure(COORD_SC, SC_VECTORS, 6, [](const Vector3& v1, const Vector3& v2){
        const double dot = v1.dot(v2);
        if(std::abs(dot) < EPSILON){
            return true;
        }
        return (v1 + v2).squaredLength() < EPSILON;
    }, [](int){ return 0; });

    initializeLatticeStructure(LATTICE_SC, SC_VECTORS, 6, &_coordinationStructures[COORD_SC]);
    _latticeStructures[LATTICE_SC].primitiveCell.column(0) = SC_PRIMITIVE_CELL[0];
    _latticeStructures[LATTICE_SC].primitiveCell.column(1) = SC_PRIMITIVE_CELL[1];
    _latticeStructures[LATTICE_SC].primitiveCell.column(2) = SC_PRIMITIVE_CELL[2];

    _latticeStructures[LATTICE_SC].permutations.clear();
    for(const auto& rotation : AnalysisSymmetryUtils::cubicSymmetryRotations()){
        SymmetryPermutation permutation;
        permutation.transformation = rotation;
        for(int i = 0; i < 6; ++i){
            const Vector3 vector = rotation * SC_VECTORS[i];
            for(int j = 0; j < 6; ++j){
                if(vector.equals(SC_VECTORS[j])){
                    permutation.permutation[i] = j;
                    break;
                }
            }
        }

        _latticeStructures[LATTICE_SC].permutations.push_back(std::move(permutation));
    }
}

void CoordinationStructures::initializeOther(){
    _coordinationStructures[COORD_OTHER].numNeighbors = 0;
    _latticeStructures[LATTICE_OTHER].coordStructure = &_coordinationStructures[COORD_OTHER];
    _latticeStructures[LATTICE_OTHER].primitiveCell.setZero();
    _latticeStructures[LATTICE_OTHER].primitiveCellInverse.setZero();
    _latticeStructures[LATTICE_OTHER].maxNeighbors = 0;
}

template <typename BondPredicate, typename SignatureFunction>
void CoordinationStructures::initializeCoordinationStructure(
    int coordType,
    const Vector3* vectors,
    int numNeighbors,
    BondPredicate bondPredicate,
    SignatureFunction signatureFunction
){
    _coordinationStructures[coordType].numNeighbors = numNeighbors;
    for(int ni1 = 0; ni1 < numNeighbors; ni1++){
        _coordinationStructures[coordType].neighborArray.setNeighborBond(ni1, ni1, false);
        for(int ni2 = ni1 + 1; ni2 < numNeighbors; ni2++){
            bool bonded = bondPredicate(vectors[ni1], vectors[ni2]);
            _coordinationStructures[coordType].neighborArray.setNeighborBond(ni1, ni2, bonded);
        }
        _coordinationStructures[coordType].cnaSignatures[ni1] = signatureFunction(ni1);
    }

    _coordinationStructures[coordType].latticeVectors.assign(vectors, vectors + numNeighbors);
}

void CoordinationStructures::initializeLatticeStructure(
    int latticeType,
    const Vector3* vectors,
    int totalVectors,
    CoordinationStructure* coordStruct
){
    _latticeStructures[latticeType].latticeVectors.assign(vectors, vectors + totalVectors);
    _latticeStructures[latticeType].coordStructure = coordStruct;
    _latticeStructures[latticeType].maxNeighbors = coordStruct->numNeighbors;
}

void CoordinationStructures::initializeDiamondStructure(
    int coordType,
    int latticeType,
    const Vector3* vectors,
    int numNeighbors,
    int totalVectors
){
    _coordinationStructures[coordType].numNeighbors = numNeighbors;
    for(int ni1 = 0; ni1 < numNeighbors; ++ni1){
        _coordinationStructures[coordType].neighborArray.setNeighborBond(ni1, ni1, false);
        double cutoff = (ni1 < 4) ? (std::sqrt(3.0) * 0.25 + std::sqrt(0.5)) / 2.0 : (1.0 + std::sqrt(0.5)) / 2.0;

        for(int ni2 = 0; ni2 < 4; ++ni2){
            if(ni1 < 4 && ni2 < 4){
                _coordinationStructures[coordType].neighborArray.setNeighborBond(ni1, ni2, false);
            }
        }

        for(int ni2 = std::max(ni1 + 1, 4); ni2 < numNeighbors; ++ni2){
            bool bonded = (vectors[ni1] - vectors[ni2]).length() < cutoff;
            _coordinationStructures[coordType].neighborArray.setNeighborBond(ni1, ni2, bonded);
        }

        if(coordType == COORD_HEX_DIAMOND){
            _coordinationStructures[coordType].cnaSignatures[ni1] = (ni1 < 4) ? 0 : ((vectors[ni1].z() == 0) ? 2 : 1);
        }else{
            _coordinationStructures[coordType].cnaSignatures[ni1] = (ni1 < 4) ? 0 : 1;
        }
    }
    _coordinationStructures[coordType].latticeVectors.assign(vectors, vectors + numNeighbors);
    initializeLatticeStructure(latticeType, vectors, totalVectors, &_coordinationStructures[coordType]);
}

void CoordinationStructures::initializeCommonNeighbors(){
    for(auto coordStruct = std::begin(_coordinationStructures); coordStruct != std::end(_coordinationStructures); ++coordStruct){
        for(int neighborIndex = 0; neighborIndex < coordStruct->numNeighbors; neighborIndex++){
            findCommonNeighborsForBond(*coordStruct, neighborIndex);
        }
    }
}

void CoordinationStructures::findCommonNeighborsForBond(CoordinationStructure& coordStruct, int neighborIndex){
    Matrix3 tm;
    tm.column(0) = coordStruct.latticeVectors[neighborIndex];
    bool found = false;

    coordStruct.commonNeighbors[neighborIndex][0] = -1;
    coordStruct.commonNeighbors[neighborIndex][1] = -1;

    if(coordStruct.numNeighbors == 6){
        for(int i1 = 0; i1 < 6 && !found; i1++){
            if(i1 == neighborIndex || i1 == (neighborIndex ^ 1)){
                continue;
            }
            tm.column(1) = coordStruct.latticeVectors[i1];

            for(int i2 = i1 + 1; i2 < 6; i2++){
                if(i2 == neighborIndex || i2 == (neighborIndex ^ 1) || i2 == (i1 ^ 1)){
                    continue;
                }
                tm.column(2) = coordStruct.latticeVectors[i2];

                if(std::abs(tm.determinant()) > EPSILON){
                    coordStruct.commonNeighbors[neighborIndex][0] = i1;
                    coordStruct.commonNeighbors[neighborIndex][1] = i2;
                    found = true;
                    break;
                }
            }
        }
        if(found){
            return;
        }
    }

    for(int i1 = 0; i1 < coordStruct.numNeighbors && !found; i1++){
        if(!coordStruct.neighborArray.neighborBond(neighborIndex, i1)){
            continue;
        }
        tm.column(1) = coordStruct.latticeVectors[i1];

        for(int i2 = i1 + 1; i2 < coordStruct.numNeighbors; i2++){
            if(!coordStruct.neighborArray.neighborBond(neighborIndex, i2)){
                continue;
            }
            tm.column(2) = coordStruct.latticeVectors[i2];

            if(std::abs(tm.determinant()) > EPSILON){
                coordStruct.commonNeighbors[neighborIndex][0] = i1;
                coordStruct.commonNeighbors[neighborIndex][1] = i2;
                found = true;
                break;
            }
        }
    }
}

void CoordinationStructures::initializeSymmetryInformation(){
    for(auto latticeStruct = std::begin(_latticeStructures); latticeStruct != std::end(_latticeStructures); ++latticeStruct){
        if(latticeStruct->latticeVectors.empty()){
            continue;
        }

        latticeStruct->primitiveCellInverse = latticeStruct->primitiveCell.inverse();
        if(latticeStruct->permutations.empty()){
            generateSymmetryPermutations(*latticeStruct);
        }
        calculateSymmetryProducts(*latticeStruct);
    }
}

void CoordinationStructures::generateSymmetryPermutations(LatticeStructure& latticeStruct){
    const CoordinationStructure& coordStruct = *latticeStruct.coordStructure;
    AnalysisSymmetryUtils::generateSymmetryPermutations(
        coordStruct.latticeVectors,
        coordStruct.numNeighbors,
        latticeStruct.latticeVectors,
        latticeStruct.permutations
    );
}

void CoordinationStructures::findNonCoplanarVectors(const CoordinationStructure& coordStruct, int nindices[3], Matrix3& tm1){
    AnalysisSymmetryUtils::findNonCoplanarVectors(coordStruct.latticeVectors, coordStruct.numNeighbors, nindices, tm1);
}

void CoordinationStructures::calculateSymmetryProducts(LatticeStructure& latticeStruct){
    AnalysisSymmetryUtils::calculateSymmetryProducts(latticeStruct.permutations);
}

void CoordinationStructures::initializeStructures(){
    initializeOther();
    initializeFCC();
    initializeHCP();
    initializeBCC();
    initializeSC();
    initializeCubicDiamond();
    initializeHexagonalDiamond();
    initializeCommonNeighbors();
    initializeSymmetryInformation();
}

}

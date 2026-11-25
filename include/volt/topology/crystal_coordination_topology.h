#pragma once

#include <volt/core/volt.h>
#include <volt/topology/crystal_symmetry_topology.h>
#include <volt/topology/crystal_coordination_pattern.h>
#include <volt/analysis/nearest_neighbor_finder.h>
#include <volt/core/particle_property.h>
#include <volt/topology/native_lattice_vectors.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/cna_local_structure.h>

namespace Volt{
    
class CoordinationStructures{
public:
    CoordinationStructures(ParticleProperty* structureTypes, LatticeStructureType inputCrystalType, bool identifyPlanarDefects, const SimulationCell& simCell);

    double determineLocalStructure(
        const NearestNeighborFinder& neighList,
        int particleIndex,
        int* outNeighborCount,
        int* outOrderedNeighborIndices = nullptr
    ) const;
    
    static void initializeStructures();

    void postProcessDiamondNeighbors(
        AnalysisContext& context,
        const NearestNeighborFinder& neighList
    ) const;

    static const LatticeStructureType getLatticeIdx(int s){
        switch(s){
            case StructureType::SC:  return LATTICE_SC;
            case StructureType::FCC: return LATTICE_FCC;
            case StructureType::HCP: return LATTICE_HCP;
            case StructureType::BCC: return LATTICE_BCC;
            case StructureType::CUBIC_DIAMOND:
            case StructureType::CUBIC_DIAMOND_FIRST_NEIGH:
            case StructureType::CUBIC_DIAMOND_SECOND_NEIGH:
                return LATTICE_CUBIC_DIAMOND;
            case StructureType::HEX_DIAMOND:
            case StructureType::HEX_DIAMOND_FIRST_NEIGH:
            case StructureType::HEX_DIAMOND_SECOND_NEIGH:
                return LATTICE_HEX_DIAMOND;
            case StructureType::ICO:
            case StructureType::GRAPHENE:
                return LATTICE_OTHER;
            default:
                if(s != StructureType::OTHER){
                    spdlog::warn("getLatticeIdx: unknown {}", s);
                }
                return LATTICE_OTHER;
        }
    }

    static const CoordinationStructureType getCoordIdx(int s){
        switch(s){
            case StructureType::SC:  return COORD_SC;
            case StructureType::FCC: return COORD_FCC;
            case StructureType::HCP: return COORD_HCP;
            case StructureType::BCC: return COORD_BCC;
            case StructureType::CUBIC_DIAMOND:
            case StructureType::CUBIC_DIAMOND_FIRST_NEIGH:
            case StructureType::CUBIC_DIAMOND_SECOND_NEIGH:
                return COORD_CUBIC_DIAMOND;
            case StructureType::HEX_DIAMOND:
            case StructureType::HEX_DIAMOND_FIRST_NEIGH:
            case StructureType::HEX_DIAMOND_SECOND_NEIGH:
                return COORD_HEX_DIAMOND;
            case StructureType::ICO:
            case StructureType::GRAPHENE:
                return COORD_OTHER;
            default:
                if(s != StructureType::OTHER){
                    spdlog::warn("getCoordIdx: unknown {}", s);
                }
                return COORD_OTHER;
        }
    }

    static const CoordinationStructure& getCoordStruct(int structureType){
        return _coordinationStructures[static_cast<std::size_t>(CoordinationStructures::getCoordIdx(structureType))];
    }

    static const LatticeStructure& getLatticeStruct(int structureType){
        return _latticeStructures[static_cast<std::size_t>(CoordinationStructures::getLatticeIdx(structureType))];
    }

    static const LatticeStructure& getLatticeStructByLatticeType(int latticeType){
        return _latticeStructures[static_cast<std::size_t>(latticeType)];
    }

    static CoordinationStructure _coordinationStructures[NUM_COORD_TYPES];
    static LatticeStructure _latticeStructures[NUM_LATTICE_TYPES];
    const SimulationCell& cell() const{
        return _simCell;
    }

    static void findNonCoplanarVectors(const CoordinationStructure& coordStruct, int nindices[3], Matrix3& tm1);

private:
    static void initializeDiamondStructure(
        int coordType,
        int latticeType,
        const Vector3* vectors,
        int numNeighbors,
        int totalVectors
    );

    static void initializeLatticeStructure(
        int latticeType,
        const Vector3* vectors,
        int totalVectors,
        CoordinationStructure* coordStruct
    );

    template <typename BondPredicate, typename SignatureFunction>
    static void initializeCoordinationStructure(
        int coordType,
        const Vector3* vectors,
        int numNeighbors,
        BondPredicate bondPredicate,
        SignatureFunction signatureFunction
    );

    static void initializeFCC();
    static void initializeSC();
    static void initializeHCP();
    static void initializeBCC();
    static void initializeCubicDiamond();
    static void initializeHexagonalDiamond();
    static void initializeOther();

    static void calculateSymmetryProducts(LatticeStructure& latticeStruct);
    static void generateSymmetryPermutations(LatticeStructure& latticeStruct);
    static void initializeSymmetryInformation();
    static void findCommonNeighborsForBond(CoordinationStructure& coordStruct, int neighborIndex);
    static void initializeCommonNeighbors();

    static bool latticeSupportsOrientationInversion(
        const LatticeStructure& lattice,
        int coordinationNumber,
        std::array<int, MAX_NEIGHBORS>& outInversePermutation
    );
    static bool computeMatchedOrientation(
        const LatticeStructure& lattice,
        const Vector3* neighborVectors,
        const int* neighborMapping,
        int coordinationNumber,
        Matrix3& outRawOrientation,
        Matrix3& outOrientation
    );
    static bool canonicalizeNeighborMapping(
        const LatticeStructure& lattice,
        const Vector3* neighborVectors,
        int coordinationNumber,
        int* neighborMapping
    );
        
    const SimulationCell& _simCell;
    ParticleProperty* _structureTypes;
    LatticeStructureType _inputCrystalType;
    bool _identifyPlanarDefects;
};

}

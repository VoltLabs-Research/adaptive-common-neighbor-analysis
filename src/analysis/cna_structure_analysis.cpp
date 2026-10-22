#include <volt/analysis/cna_structure_analysis.h>
#include <volt/analysis/crystal_symmetry_utils.h>
#include <volt/core/progress.h>
#include <volt/topology/crystal_coordination_topology.h>
#include <volt/topology/crystal_coordination_topology_init.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <algorithm>
#include <vector>

namespace Volt {

namespace CnaStructureAnalysisDetail{

class CnaCrystalInfoProvider final : public StructureAnalysisCrystalInfo{
public:
    std::string_view topologyName(int structureType) const override{
        switch(structureType){
            case LATTICE_SC:
                return "sc";
            case LATTICE_FCC:
                return "fcc";
            case LATTICE_HCP:
                return "hcp";
            case LATTICE_BCC:
                return "bcc";
            case LATTICE_CUBIC_DIAMOND:
                return "cubic_diamond";
            case LATTICE_HEX_DIAMOND:
                return "hex_diamond";
            default:
                return {};
        }
    }

    int findClosestSymmetryPermutation(int structureType, const Matrix3& rotation) const override{
        ensureCoordinationStructuresInitialized();

        const LatticeStructure& lattice = CoordinationStructures::getLatticeStruct(structureType);
        return AnalysisSymmetryUtils::findClosestSymmetryPermutation(lattice.permutations, rotation);
    }

    int coordinationNumber(int structureType) const override{
        ensureCoordinationStructuresInitialized();
        return CoordinationStructures::getCoordStruct(structureType).numNeighbors;
    }

    int commonNeighborIndex(int structureType, int neighborIndex, int commonNeighborSlot) const override{
        ensureCoordinationStructuresInitialized();
        return CoordinationStructures::getCoordStruct(structureType)
            .commonNeighbors[neighborIndex][commonNeighborSlot];
    }

    int symmetryPermutationCount(int structureType) const override{
        ensureCoordinationStructuresInitialized();
        return static_cast<int>(CoordinationStructures::getLatticeStruct(structureType).permutations.size());
    }

    int symmetryPermutationEntry(int structureType, int symmetryIndex, int neighborIndex) const override{
        ensureCoordinationStructuresInitialized();
        return CoordinationStructures::getLatticeStruct(structureType)
            .permutations[symmetryIndex].permutation[neighborIndex];
    }

    const Matrix3& symmetryTransformation(int structureType, int symmetryIndex) const override{
        ensureCoordinationStructuresInitialized();
        return CoordinationStructures::getLatticeStruct(structureType)
            .permutations[symmetryIndex].transformation;
    }

    int symmetryInverseProduct(int structureType, int symmetryIndex, int transformationIndex) const override{
        ensureCoordinationStructuresInitialized();
        return CoordinationStructures::getLatticeStruct(structureType)
            .permutations[symmetryIndex].inverseProduct[transformationIndex];
    }

    const Vector3& latticeVector(int structureType, int latticeVectorIndex) const override{
        ensureCoordinationStructuresInitialized();
        return CoordinationStructures::getLatticeStruct(structureType).latticeVectors[latticeVectorIndex];
    }
};

std::shared_ptr<const StructureAnalysisCrystalInfo> cnaCrystalInfoProvider(){
    static const auto provider = std::make_shared<CnaCrystalInfoProvider>();
    return provider;
}

}

using namespace CnaStructureAnalysisDetail;

void identifyStructuresCNA(StructureAnalysis& analysis){
    StructureContext& context = analysis.context();
    ensureCoordinationStructuresInitialized();
    analysis.setCrystalInfoProvider(cnaCrystalInfoProvider());

    const Stopwatch neighborTimer;
    const int maxNeighborListSize = MAX_NEIGHBORS;
    NearestNeighborFinder neighFinder(maxNeighborListSize);
    if(!neighFinder.prepare(context.positions, context.simCell)){
        throw std::runtime_error("Error in neighFinder.preapre(...)");
    }
    spdlog::info("  Nearest-neighbor search ready in {}", neighborTimer.elapsed());

    CoordinationStructures coordinationStructures(
        context.structureTypes,
        context.inputCrystalType,
        true,
        context.simCell
    );

    const size_t N = context.atomCount();
    std::vector<int> localCounts(N, 0);
    std::vector<std::array<int, MAX_NEIGHBORS>> orderedNeighborIndices(N);
    ProgressReporter progress("Classifying atoms", N);
    context.maximumNeighborDistance = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, N),
        0.0,
        [&neighFinder, &coordinationStructures, &localCounts, &orderedNeighborIndices, &progress](const tbb::blocked_range<size_t>& range, double maxDistSoFar) -> double {
            for(size_t index = range.begin(); index != range.end(); ++index){
                int count = 0;
                auto& neighbors = orderedNeighborIndices[index];
                neighbors.fill(-1);
                double localMaxDistance = coordinationStructures.determineLocalStructure(
                    neighFinder,
                    index,
                    &count,
                    neighbors.data()
                );
                localCounts[index] = count;
                if(localMaxDistance > maxDistSoFar){
                    maxDistSoFar = localMaxDistance;
                }
            }
            progress.advance(range.size());
            return maxDistSoFar;
        },
        [](double a, double b) -> double {
            return std::max(a, b);
        }
    );

    auto* offsets = context.neighborOffsets->dataInt();
    offsets[0] = 0;
    for(size_t i = 0; i < N; ++i){
        offsets[i + 1] = offsets[i] + localCounts[i];
    }
    const size_t totalNeighbors = static_cast<size_t>(offsets[N]);
    context.neighborIndices = std::make_shared<ParticleProperty>(
        totalNeighbors, DataType::Int, 1, 0, false);

    auto* indices = context.neighborIndices->dataInt();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, N), [&](const auto& range){
        for(size_t index = range.begin(); index != range.end(); ++index){
            const int count = localCounts[index];
            if(count <= 0){
                continue;
            }
            const int start = offsets[index];
            for(int j = 0; j < count; ++j){
                indices[start + j] = orderedNeighborIndices[index][j];
            }
            context.neighborCounts->setInt(index, count);
        }
    });
}

}

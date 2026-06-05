#include <volt/analysis/cna_analysis.h>
#include <volt/core/neighbor_ordering.h>
#include <volt/topology/crystal_coordination_topology_init.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <cstdint>

namespace Volt{

CommonNeighborAnalysisEngine::CommonNeighborAnalysisEngine(AnalysisContext& context, bool identifyPlanarDefects)
    : _context(context), _identifyPlanarDefects(identifyPlanarDefects){}

void CommonNeighborAnalysisEngine::perform(){
    identifyStructures();
    calculateStructureStatistics();
}

double CommonNeighborAnalysisEngine::determineLocalStructure(
    const NearestNeighborFinder& neighList,
    int particleIndex
) const{
    CnaLocalStructureUtils::LocalStructureMatch match;
    if(!CnaLocalStructureUtils::determineLocalStructure(
        _context.inputCrystalType,
        _identifyPlanarDefects,
        neighList,
        particleIndex,
        CoordinationStructures::_coordinationStructures,
        match
    )){
        return 0.0;
    }

    _context.structureTypes->setInt(particleIndex, static_cast<int>(match.structureType));
    storeNeighborOrdering(
        match.neighborMapping,
        match.coordinationNumber,
        match.structureType,
        particleIndex
    );

    return match.localCutoff;
}

void CommonNeighborAnalysisEngine::storeNeighborOrdering(
    const std::array<int, MAX_NEIGHBORS>& neighborMapping,
    int coordinationCount,
    StructureType atomStructure,
    int particleIndex
) const{
    if(!_context.correspondences || coordinationCount <= 0){
        return;
    }

    if(!describeNeighborOrdering(atomStructure).isSingleShell()){
        return;
    }

    NeighborOrderingStorage encodedOrdering{};
    encodedOrdering[0] = 0;
    for(int i = 0; i < coordinationCount; ++i){
        encodedOrdering[static_cast<size_t>(i + 1)] = static_cast<int8_t>(neighborMapping[static_cast<size_t>(i)] + 1);
    }

    std::uint64_t orderingCode = 0;
    if(!encodeNeighborOrdering(atomStructure, coordinationCount, encodedOrdering, orderingCode)){
        return;
    }

    _context.correspondences->setInt64(static_cast<size_t>(particleIndex), static_cast<std::int64_t>(orderingCode));
}

void CommonNeighborAnalysisEngine::identifyStructures(){
    ensureCoordinationStructuresInitialized();

    std::fill(
        _context.structureTypes->dataInt(),
        _context.structureTypes->dataInt() + _context.structureTypes->size(),
        static_cast<int>(StructureType::OTHER)
    );

    NearestNeighborFinder neighFinder(MAX_NEIGHBORS);
    if(!neighFinder.prepare(_context.positions, _context.simCell)){
        throw std::runtime_error("Error in neighFinder.prepare(...)");
    }

    const size_t N = _context.atomCount();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, N), [&](const tbb::blocked_range<size_t>& r){
        for(size_t index = r.begin(); index != r.end(); ++index){
            determineLocalStructure(neighFinder, static_cast<int>(index));
        }
    });
}

void CommonNeighborAnalysisEngine::calculateStructureStatistics() const{
}

std::string CommonNeighborAnalysisEngine::getStructureTypeName(int structureType) const{
    return std::string(structureTypeName(structureType));
}

}


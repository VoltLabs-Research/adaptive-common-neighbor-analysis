#include <volt/analysis/common_neighbor_analysis_log.h>

#include <volt/analysis/structure_analysis.h>
#include <volt/core/type_counts.h>

#include <spdlog/spdlog.h>

namespace Volt{

CommonNeighborAnalysisLog::CommonNeighborAnalysisLog(LatticeStructureType inputLattice)
    : _assignableTypes(assignableStructureTypes(inputLattice)){}

std::vector<int> CommonNeighborAnalysisLog::assignableStructureTypes(LatticeStructureType lattice){
    switch(lattice){
        case LATTICE_FCC:
        case LATTICE_HCP:
            return {StructureType::FCC, StructureType::HCP};
        case LATTICE_BCC:
            return {StructureType::BCC};
        case LATTICE_CUBIC_DIAMOND:
        case LATTICE_HEX_DIAMOND:
            return {StructureType::CUBIC_DIAMOND, StructureType::HEX_DIAMOND};
        default:
            return {};
    }
}

std::string CommonNeighborAnalysisLog::signatureNames(const std::vector<int>& structureTypes){
    std::string names;
    for(const int structureType : structureTypes){
        if(!names.empty()) names += '/';
        names += structureTypeName(structureType);
    }
    return names;
}

void CommonNeighborAnalysisLog::dump(const StructureAnalysis& classified){
    const StructureContext& context = classified.context();
    spdlog::info(
        "Classified {} atoms in {} (maximum adaptive cutoff {:.3f})",
        context.atomCount(), phaseElapsed(), context.maximumNeighborDistance
    );
}

void CommonNeighborAnalysisLog::dump(const std::map<int, int>& structureCounts){
    std::map<int, int> counts = structureCounts;
    for(const int structureType : _assignableTypes){
        counts.try_emplace(structureType, 0);
    }
    TypeCounts line{"Structure types"};
    for(const auto& [structureType, count] : counts){
        line.atoms += static_cast<std::size_t>(count);
        if(structureType != StructureType::OTHER){
            line.counts.emplace_back(structureTypeName(structureType), count);
        }
    }
    const int otherAtoms = counts[StructureType::OTHER];
    line.counts.emplace_back("Other", otherAtoms);
    dump(line);

    if(line.atoms > 0 && static_cast<std::size_t>(otherAtoms) == line.atoms){
        spdlog::warn(
            "No atom matched the {} signatures; check that crystal_structure fits the sample",
            signatureNames(_assignableTypes)
        );
    }
}

}

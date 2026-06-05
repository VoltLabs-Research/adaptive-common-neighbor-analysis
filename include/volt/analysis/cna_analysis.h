#pragma once

#include <volt/core/volt.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/cna_local_structure.h>
#include <volt/core/lammps_parser.h>
#include <array>

namespace Volt{

class CommonNeighborAnalysisEngine{
public:
    CommonNeighborAnalysisEngine(AnalysisContext& context, bool identifyPlanarDefects);

    void perform();

    std::string getStructureTypeName(int structureType) const;

private:
    void identifyStructures();
    void calculateStructureStatistics() const;

    double determineLocalStructure(
        const NearestNeighborFinder& neighList,
        int particleIndex
    ) const;

    void storeNeighborOrdering(
        const std::array<int, MAX_NEIGHBORS>& neighborMapping,
        int coordinationCount,
        StructureType atomStructure,
        int particleIndex
    ) const;

    AnalysisContext& _context;
    bool _identifyPlanarDefects = false;
};

}

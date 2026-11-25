#pragma once

#include <volt/analysis/structure_identification_log.h>
#include <volt/structures/crystal_structure_types.h>

#include <map>
#include <string>
#include <vector>

namespace Volt{

class StructureAnalysis;

// The execution log of one adaptive CNA run, on top of the lines every structure identification prints.
class CommonNeighborAnalysisLog : public StructureIdentificationLog{
public:
    explicit CommonNeighborAnalysisLog(LatticeStructureType inputLattice);

    using StructureIdentificationLog::dump;

    void dump(const StructureAnalysis& classified);
    // The atoms of each structure type, keyed by type.
    void dump(const std::map<int, int>& structureCounts);

private:
    // CNA always identifies planar defects, so an FCC or HCP input lattice yields both FCC and HCP atoms, and
    // either diamond lattice yields both diamond types.
    static std::vector<int> assignableStructureTypes(LatticeStructureType lattice);
    static std::string signatureNames(const std::vector<int>& structureTypes);

    // Listed in the structure counts even when no atom has them.
    const std::vector<int> _assignableTypes;
};

}

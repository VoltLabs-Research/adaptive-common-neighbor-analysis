#include <volt/plugin/plugin_entry.h>
#include <volt/analysis/cna_service.h>
#include <volt/structures/crystal_structure_types.h>

using namespace Volt;
using namespace Volt::Plugin;
using S = CommonNeighborAnalysisService;

static const std::vector<OptionBinding<S>> bindings = {
    opt<S>({"--crystal_structure", "string", "FCC|BCC|HCP|CUBIC_DIAMOND|HEX_DIAMOND", "FCC"},
        [](S& s, const OptsMap& opts) {
            LatticeStructureType t = LATTICE_FCC;
            parseLatticeStructureType(CLI::getString(opts, "--crystal_structure", "FCC"), t);
            s.setInputCrystalStructure(t);
        }),
    opt("--dissolve_small_clusters", "Mark small clusters as OTHER", false, &S::setDissolveSmallClusters),
};

VOLT_SERVICE_PLUGIN("volt-common-neighbor-analysis", "Common Neighbor Analysis", S, bindings)

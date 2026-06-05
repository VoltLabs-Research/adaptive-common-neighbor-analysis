#include <volt/plugin/plugin_entry.h>
#include <volt/analysis/cna_service.h>

using namespace Volt;
using namespace Volt::Plugin;
using S = CommonNeighborAnalysisService;

static const std::vector<OptionBinding<S>> bindings = {
    optLattice("--crystal_structure", "FCC|BCC|HCP|CUBIC_DIAMOND|HEX_DIAMOND", "FCC", &S::setInputCrystalStructure),
    opt("--dissolve_small_clusters", "Mark small clusters as OTHER", false, &S::setDissolveSmallClusters),
};

VOLT_SERVICE_PLUGIN("volt-common-neighbor-analysis", "Common Neighbor Analysis", S, bindings)

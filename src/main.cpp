#include <volt/cli/common.h>
#include <volt/analysis/cna_service.h>
#include <volt/structures/crystal_structure_types.h>

#include <string>

using namespace Volt;
using namespace Volt::CLI;

namespace Volt::CnaCliDetail {

void showUsage(const std::string& name){
    printUsageHeader(name, "Volt - Common Neighbor Analysis");
    std::cerr
        << "  --crystal_structure <type>     Crystal structure. (FCC|BCC|HCP|CUBIC_DIAMOND|HEX_DIAMOND) [default: FCC]\n"
        << "  --dissolve_small_clusters       Mark small clusters as OTHER after building clusters.\n";
    printHelpOption();
}

}

int main(int argc, char* argv[]){
    std::string filename, outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);
    if(const int startupStatus = handleHelpOrMissingInput(argc, argv, opts, filename, Volt::CnaCliDetail::showUsage);
       startupStatus >= 0){
        return startupStatus;
    }

    initLogging("volt-common-neighbor-analysis");

    LammpsParser::Frame frame;
    if(!parseFrame(filename, frame)) return 1;

    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);

    CommonNeighborAnalysisService analyzer;
    LatticeStructureType crystal_structure = LATTICE_FCC;
    const std::string crystal_structureOption = getString(opts, "--crystal_structure", "FCC");
    if(!parseLatticeStructureType(crystal_structureOption, crystal_structure)){
        spdlog::warn("Unknown crystal structure '{}', defaulting to FCC.", crystal_structureOption);
        crystal_structure = LATTICE_FCC;
    }
    analyzer.setInputCrystalStructure(crystal_structure);
    analyzer.setDissolveSmallClusters(hasOption(opts, "--dissolve_small_clusters"));

    spdlog::info("Starting common neighbor analysis...");
    json result = analyzer.compute(frame, outputBase, filename);
    if(result.value("is_failed", false)){
        spdlog::error("Analysis failed: {}", result.value("error", "Unknown error"));
        return 1;
    }

    spdlog::info("Common neighbor analysis completed.");
    return 0;
}

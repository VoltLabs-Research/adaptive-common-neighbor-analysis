#include <volt/analysis/cna_service.h>

#include <volt/analysis/reconstructed_analysis_pipeline.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/cluster_graph_builder.h>
#include <volt/analysis/cluster_graph_io.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/structure_identification_export.h>
#include <volt/analysis/cna_cluster_input_adapter.h>
#include <volt/analysis/cna_local_structure.h>
#include <volt/analysis/cna_structure_analysis.h>
#include <volt/core/analysis_result.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/particle_property.h>
#include <volt/core/progress.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/utilities/json_utils.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <map>
#include <utility>
#include <vector>

namespace Volt{

using namespace Volt::Particles;

namespace{

// identifyStructuresCNA always identifies planar defects, so an FCC or HCP input lattice yields
// both FCC and HCP atoms, and either diamond lattice yields both diamond types.
std::vector<int> assignableStructureTypes(LatticeStructureType lattice){
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

}


CommonNeighborAnalysisService::CommonNeighborAnalysisService()
    : _inputCrystalStructure(LATTICE_FCC)
    , _dissolveSmallClusters(false){}

void CommonNeighborAnalysisService::setInputCrystalStructure(LatticeStructureType structureType){
    _inputCrystalStructure = structureType;
}

void CommonNeighborAnalysisService::setDissolveSmallClusters(bool dissolveSmallClusters){
    _dissolveSmallClusters = dissolveSmallClusters;
}

json CommonNeighborAnalysisService::compute(
    const LammpsParser::Frame& frame,
    const std::string& outputBase,
    const std::string& inputDumpPath
){
    const std::string annotatedDumpPath = outputBase.empty()
        ? inputDumpPath
        : outputBase + "_annotated.dump";

    if(_inputCrystalStructure == LATTICE_SC){
        return AnalysisResult::failure("CNA does not support SC. Use PTM for this crystal.");
    }

    std::string frameError;
    auto session = AnalysisPipelineUtils::prepareAnalysisSession(
        frame,
        _inputCrystalStructure,
        &frameError
    );
    if(!session){
        return AnalysisResult::failure(frameError);
    }
    AnalysisContext& context = session->context;

    try{
        StructureAnalysis analysis(context);
        const std::vector<int> assignableTypes = assignableStructureTypes(_inputCrystalStructure);
        std::string signatures;
        for(const int structureType : assignableTypes){
            if(!signatures.empty()) signatures += '/';
            signatures += structureTypeName(structureType);
        }
        spdlog::info(
            "Classifying {} atoms with adaptive CNA ({} signatures, {}-neighbor shell)...",
            context.atomCount(),
            signatures,
            CnaLocalStructureUtils::coordinationNumberFor(_inputCrystalStructure)
        );
        Stopwatch timer;
        identifyStructuresCNA(analysis);
        spdlog::info(
            "Classified {} atoms in {} (maximum adaptive cutoff {:.3f})",
            context.atomCount(),
            timer.elapsed(),
            context.maximumNeighborDistance
        );

        CNAClusterInputAdapter clusterInputAdapter;
        clusterInputAdapter.prepare(analysis, context);
        ClusterBuilder clusterBuilder(analysis, context);
        clusterBuilder.build(_dissolveSmallClusters);

        std::map<int,int> structCounts;
        for(int i = 0; i < frame.natoms; ++i)
            structCounts[context.structureTypes->getInt(i)]++;

        json result;
        result["main_listing"] = {
            {"total_atoms", frame.natoms},
            {"structure_count", static_cast<int>(structCounts.size())}
        };
        result["sub_listings"] = json::object();
        result["per-atom-properties"] = json::array();

        if(!outputBase.empty()){
            const std::string parquetPath = outputBase + "_cna_analysis.parquet";
            JsonUtils::writeJsonToParquet(result, parquetPath, false);

            const std::string atomsPath = outputBase + "_atoms.parquet";
            StructureIdentificationExport::streamStructureIdentificationToParquet(
                atomsPath, frame, analysis
            );
        }

        if(!AnalysisPipelineUtils::appendClusterOutputs(
            frame,
            outputBase,
            annotatedDumpPath,
            context,
            analysis,
            result,
            &frameError
        )){
            return AnalysisResult::failure(frameError);
        }

        result["is_failed"] = false;
        return result;
    }catch(const std::exception& error){
        return AnalysisResult::failure(std::string("CNA analysis failed: ") + error.what());
    }
}

}

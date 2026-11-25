#include <volt/analysis/cna_service.h>

#include <volt/analysis/reconstructed_analysis_pipeline.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/cluster_graph_builder.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/structure_identification_export.h>
#include <volt/analysis/cna_cluster_input_adapter.h>
#include <volt/analysis/cna_structure_analysis.h>
#include <volt/analysis/common_neighbor_analysis_log.h>
#include <volt/core/analysis_result.h>
#include <volt/core/progress.h>
#include <volt/core/type_counts.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/utilities/json_utils.h>

#include <map>

namespace Volt{

using namespace Volt::Particles;

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
        CommonNeighborAnalysisLog log(_inputCrystalStructure);
        StructureAnalysis analysis(context);
        ProgressSteps steps({
            {"Classifying atoms", 1.0},
            {"Building clusters", 3.0},
            {"Writing the atoms table", 1.0},
            {"Writing cluster outputs", 4.0},
        });
        identifyStructuresCNA(analysis);
        log.dump(analysis);
        steps.next();

        CNAClusterInputAdapter clusterInputAdapter;
        clusterInputAdapter.prepare(analysis, context);
        log.dump(ClusterBuilder(analysis, context).build(_dissolveSmallClusters));

        const std::map<int, int> structCounts = countTypes(context.structureTypes->constIntRange());
        log.dump(structCounts);

        json result;
        result["main_listing"] = {
            {"total_atoms", frame.natoms},
            {"structure_count", static_cast<int>(structCounts.size())}
        };
        result["sub_listings"] = json::object();
        result["per-atom-properties"] = json::array();
        steps.next();

        if(!outputBase.empty()){
            JsonUtils::writeSummaryTable(result, outputBase, "cna_analysis");
            log.dump(StructureIdentificationExport::writeAtomsTable(outputBase, frame, analysis));
        }
        steps.next();

        const auto clusterOutputs = AnalysisPipelineUtils::appendClusterOutputs(
            frame,
            outputBase,
            annotatedDumpPath,
            context,
            analysis,
            result,
            &frameError
        );
        if(!clusterOutputs){
            return AnalysisResult::failure(frameError);
        }
        log.dump(clusterOutputs);

        result["is_failed"] = false;
        return result;
    }catch(const std::exception& error){
        return AnalysisResult::failure(std::string("CNA analysis failed: ") + error.what());
    }
}

}

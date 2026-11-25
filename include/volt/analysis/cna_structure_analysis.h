#pragma once

#include <volt/analysis/structure_analysis.h>

#include <memory>

namespace Volt{

void identifyStructuresCNA(StructureAnalysis& analysis);

namespace CnaStructureAnalysisDetail{

// The crystal information of the CNA structure types, one instance shared by every analysis.
std::shared_ptr<const StructureAnalysisCrystalInfo> cnaCrystalInfoProvider();

}

}

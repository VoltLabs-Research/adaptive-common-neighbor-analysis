cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(common-neighbor-analysis VERSION 1.0.0 LANGUAGES CXX)

if(NOT TARGET coretoolkit::coretoolkit)
    find_package(coretoolkit REQUIRED)
endif()
if(NOT TARGET yaml-cpp::yaml-cpp)
    find_package(yaml-cpp REQUIRED)
endif()
if(NOT TARGET structure-identification::structure-identification)
    find_package(structure-identification REQUIRED)
endif()

volt_add_plugin(common-neighbor-analysis
    SOURCES
        src/analysis/cna_service.cpp
        src/analysis/cna_classifier.cpp
        src/analysis/cna_cluster_input_adapter.cpp
        src/analysis/cna_structure_analysis.cpp
        src/topology/crystal_coordination_topology.cpp
    MAIN_SOURCE src/plugin.cpp
    DEPENDENCIES
        structure-identification::structure-identification
        yaml-cpp::yaml-cpp
)

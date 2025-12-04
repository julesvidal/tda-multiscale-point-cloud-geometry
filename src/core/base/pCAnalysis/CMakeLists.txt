set (CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 17)

ttk_add_base_library(pCAnalysis
  SOURCES
    PCAnalysis.cpp
    Sampling.cpp
    RIMLS.cpp
  HEADERS
    PCAnalysisUtils.h
    RIMLSWeightFunc.h
    RIMLS.h
    PCAnalysis.h
    Sampling.h
  DEPENDS
    triangulation
)

add_subdirectory(ponca)
set_target_properties(pCAnalysis PROPERTIES LINKER_LANGUAGE CXX)

target_compile_definitions(pCAnalysis PRIVATE TTK_ENABLE_EIGEN)
target_link_libraries(pCAnalysis PUBLIC Eigen3::Eigen Ponca::Fitting)


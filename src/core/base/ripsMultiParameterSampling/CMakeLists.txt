ttk_add_base_library(ripsMultiParameterSampling
  SOURCES
    RipsMultiParameterSampling.cpp
  HEADERS
    RipsMultiParameterSampling.h
  DEPENDS
  ripsPersistence
  ripsTree
  pCAnalysis
    triangulation
)


if(TTK_ENABLE_EIGEN)
  target_compile_definitions(ripsMultiParameterSampling PRIVATE TTK_ENABLE_EIGEN)
  target_link_libraries(ripsMultiParameterSampling PUBLIC Eigen3::Eigen)
endif()

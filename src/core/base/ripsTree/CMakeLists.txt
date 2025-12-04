set (CMAKE_BUILD_TYPE Release)
ttk_add_base_library(ripsTree
  SOURCES
    RipsTree.cpp
  HEADERS
    RipsTree.h
  DEPENDS
    triangulation
    persistenceDiagram
)


if(TTK_ENABLE_EIGEN)
  target_compile_definitions(ripsTree PRIVATE TTK_ENABLE_EIGEN)
  target_link_libraries(ripsTree PUBLIC Eigen3::Eigen)
endif()

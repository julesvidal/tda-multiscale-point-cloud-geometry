set (CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)

ttk_add_base_library(ripsPersistence
  SOURCES
    RipsPersistence.cpp
  HEADERS
    RipsPersistence.h
  DEPENDS
    triangulation
    persistenceDiagram
    dimensionReduction
    pCAnalysis
)

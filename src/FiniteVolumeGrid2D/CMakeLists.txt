set(HEADERS FiniteVolumeGrid2D.h
            EquidistantGrid2D.h)

set(SOURCES FiniteVolumeGrid2D.cpp
            EquidistantGrid2D.cpp)

add_library(FiniteVolumeGrid2D ${HEADERS} ${SOURCES})
target_link_libraries(FiniteVolumeGrid2D InputOutput Geometry)

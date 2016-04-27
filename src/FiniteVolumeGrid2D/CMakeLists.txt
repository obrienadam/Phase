set(HEADERS FiniteVolumeGrid2D.h
            RectilinearGrid2D.h
            ConstructGrid.h)

set(SOURCES FiniteVolumeGrid2D.cpp
            RectilinearGrid2D.cpp
            ConstructGrid.cpp)

add_subdirectory(Node)
add_subdirectory(Cell)
add_subdirectory(Face)
add_subdirectory(Link)

add_library(FiniteVolumeGrid2D ${HEADERS} ${SOURCES})
target_link_libraries(FiniteVolumeGrid2D System Geometry Cell Face Link Node)

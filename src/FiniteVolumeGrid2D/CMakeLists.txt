set(HEADERS Node.h
            Cell.h
            Face.h
            Link.h
	    Patch.h
            CellGroup.h
            FiniteVolumeGrid2D.h
            RectilinearGrid2D.h
            ConstructGrid.h
            RangeSearch.h)

set(SOURCES Node.cpp
            Cell.cpp
            Face.cpp
            Link.cpp
            Patch.cpp
            CellGroup.cpp
            FiniteVolumeGrid2D.cpp
            RectilinearGrid2D.cpp
            ConstructGrid.cpp
            RangeSearch.cpp)

add_library(FiniteVolumeGrid2D ${HEADERS} ${SOURCES})
target_link_libraries(FiniteVolumeGrid2D System Geometry)

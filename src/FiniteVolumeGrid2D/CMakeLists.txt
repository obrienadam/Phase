set(HEADERS FiniteVolumeGrid2D.h
        StructuredRectilinearGrid.h
        CgnsUnstructuredGrid.h
        ConstructGrid.h
        Node/Node.h
        Face/Face.h
        Face/Patch.h
        Cell/Cell.h
        Cell/CellGroup.h
        Cell/CellZone.h
        Link/Link.h
        Link/BoundaryFaceLink.h
        Link/InteriorFaceLink.h
        Link/DiagonalCellLink.h)

set(SOURCES FiniteVolumeGrid2D.cpp
        FiniteVolumeGrid2D.tpp
        StructuredRectilinearGrid.cpp
        CgnsUnstructuredGrid.cpp
        ConstructGrid.cpp
        Node/Node.cpp
        Face/Face.cpp
        Face/Patch.cpp
        Cell/Cell.cpp
        Cell/CellGroup.cpp
        Cell/CellZone.cpp
        Link/Link.cpp
        Link/BoundaryFaceLink.cpp
        Link/InteriorFaceLink.cpp
        Link/DiagonalCellLink.cpp)

add_library(FiniteVolumeGrid2D ${HEADERS} ${SOURCES})
target_link_libraries(FiniteVolumeGrid2D System Communicator metis Geometry)

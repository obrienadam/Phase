set(HEADERS FiniteVolumeGrid2D.h
        StructuredRectilinearGrid.h
        CgnsUnstructuredGrid.h
        ConstructGrid.h
        Node/Node.h
        Node/NodeGroup.h
        Group.h
        Face/Face.h
        Face/FaceGroup.h
        Face/Patch.h
        Cell/Cell.h
        Cell/CellGroup.h
        Cell/CellZone.h
        Link/Link.h
        Link/BoundaryFaceLink.h
        Link/InteriorFaceLink.h
        Link/DiagonalCellLink.h
        FiniteVolumeZone.h)

set(SOURCES FiniteVolumeGrid2D.cpp
        FiniteVolumeGrid2D.tpp
        StructuredRectilinearGrid.cpp
        CgnsUnstructuredGrid.cpp
        ConstructGrid.cpp
        Node/Node.cpp
        Group.tpp
        Face/Face.cpp
        Face/Patch.cpp
        Cell/Cell.cpp
        Cell/CellZone.cpp
        Link/Link.cpp
        Link/BoundaryFaceLink.cpp
        Link/InteriorFaceLink.cpp
        Link/DiagonalCellLink.cpp
        FiniteVolumeZone.cpp)

add_library(FiniteVolumeGrid2D ${HEADERS} ${SOURCES})
target_link_libraries(FiniteVolumeGrid2D metis cgns hdf5)

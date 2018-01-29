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
        Link/CellLink.h
        Link/BoundaryLink.h
        Link/InteriorLink.h
        FiniteVolumeZone.h
        BilinearInterpolator.h
        QuadraticInterpolator.h)

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
        Link/CellLink.cpp
        Link/BoundaryLink.cpp
        Link/InteriorLink.cpp
        FiniteVolumeZone.cpp
        BilinearInterpolator.cpp
        QuadraticInterpolator.cpp)

add_library(FiniteVolumeGrid2D ${HEADERS} ${SOURCES})
target_link_libraries(FiniteVolumeGrid2D Communicator Geometry System metis cgns hdf5)

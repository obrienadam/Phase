set(HEADERS FiniteVolumeGrid2D.h
        StructuredRectilinearGrid.h
        CgnsUnstructuredGrid.h
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
        BilinearInterpolator.h
        QuadraticInterpolator.h)

set(SOURCES FiniteVolumeGrid2D.cpp
        FiniteVolumeGrid2D.tpp
        StructuredRectilinearGrid.cpp
        CgnsUnstructuredGrid.cpp
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
        BilinearInterpolator.cpp
        QuadraticInterpolator.cpp FiniteVolumeGrid2DFactory.cpp FiniteVolumeGrid2DFactory.h)

add_library(FiniteVolumeGrid2D ${HEADERS} ${SOURCES})
target_link_libraries(FiniteVolumeGrid2D Communicator Geometry System metis cgns ${HDF5_LIBRARIES})

install(TARGETS
        FiniteVolumeGrid2D
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
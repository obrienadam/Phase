set(HEADERS Link.h
            BoundaryFaceLink.h
            InteriorFaceLink.h
            DiagonalCellLink.h)

set(SOURCES Link.cpp
            BoundaryFaceLink.cpp
            InteriorFaceLink.cpp
            DiagonalCellLink.cpp)

add_library(Link ${HEADERS} ${SOURCES})

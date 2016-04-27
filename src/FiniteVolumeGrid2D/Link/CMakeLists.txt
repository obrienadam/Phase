set(HEADERS Link.h
	    BoundaryLink.h
	    InteriorLink.h)

set(SOURCES Link.cpp
	    BoundaryLink.cpp
	    InteriorLink.cpp)

add_library(Link ${HEADERS} ${SOURCES})

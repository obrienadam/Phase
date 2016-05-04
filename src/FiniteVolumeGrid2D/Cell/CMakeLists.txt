set(HEADERS Cell.h
            CellGroup.h
            UniqueCellGroup.h
            CellSearch.h)

set(SOURCES Cell.cpp
            CellGroup.cpp
            UniqueCellGroup.cpp
            CellSearch.cpp)

add_library(Cell ${HEADERS} ${SOURCES})

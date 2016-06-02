set(HEADERS Cell.h
            CellGroup.h
            UniqueCellGroup.h)

set(SOURCES Cell.cpp
            CellGroup.cpp
            UniqueCellGroup.cpp)

add_library(Cell ${HEADERS} ${SOURCES})

set(HEADERS Cell.h
            CellGroup.h
            CellZone.h)

set(SOURCES Cell.cpp
            CellGroup.cpp
            CellZone.cpp)

add_library(Cell ${HEADERS} ${SOURCES})

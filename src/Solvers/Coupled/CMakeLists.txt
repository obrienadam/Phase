set(HEADERS Coupled.h
            CoupledEquation.h)

set(SOURCES Coupled.cpp
            CoupledEquation.cpp)

add_library(Coupled ${HEADERS} ${SOURCES})
target_link_libraries(Coupled Solvers)

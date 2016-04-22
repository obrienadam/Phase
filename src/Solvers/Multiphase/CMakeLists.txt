set(HEADERS Multiphase.h
            Cicsam.h
            HeightFunction.h)

set(SOURCES Multiphase.cpp
            Cicsam.cpp
            HeightFunction.cpp)

add_library(Multiphase ${HEADERS} ${SOURCES})
target_link_libraries(Multiphase Solvers)

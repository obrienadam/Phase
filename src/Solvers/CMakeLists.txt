set(HEADERS Solver.h
            Simple.h)
set(SOURCES Solver.cpp
            Simple.cpp)

add_library(Solvers ${HEADERS} ${SOURCES})
target_link_libraries(Solvers Equation)

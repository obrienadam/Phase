set(HEADERS Solver.h
            Poisson.h
            #Simple.h
            ConstructSolver.h)

set(SOURCES Solver.cpp
            Poisson.cpp
            #Simple.cpp
            ConstructSolver.cpp)

add_library(Solvers ${HEADERS} ${SOURCES})
target_link_libraries(Solvers Equation)

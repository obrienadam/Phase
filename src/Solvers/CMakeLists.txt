set(HEADERS Solver.h
            Poisson.h
            Simple.h
            Piso.h
            ConstructSolver.h)

set(SOURCES Solver.cpp
            Poisson.cpp
            Simple.cpp
            Piso.cpp
            ConstructSolver.cpp)

add_subdirectory(Multiphase)

add_library(Solvers ${HEADERS} ${SOURCES})
target_link_libraries(Solvers Equation)

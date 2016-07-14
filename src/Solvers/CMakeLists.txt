set(HEADERS Solver.h
            Poisson.h
            Simple.h
            Piso.h
            FractionalStep.h
            ConstructSolver.h)

set(SOURCES Solver.cpp
            Poisson.cpp
            Simple.cpp
            Piso.cpp
            FractionalStep.cpp
            ConstructSolver.cpp)

add_subdirectory(Multiphase)
add_subdirectory(ImmersedBoundary)
add_subdirectory(Coupled)

add_library(Solvers ${HEADERS} ${SOURCES})
target_link_libraries(Solvers Discretization)

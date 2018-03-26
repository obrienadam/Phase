set(HEADERS Solver.h
        Poisson.h
        FractionalStep.h
        FractionalStepAxisymmetric.h
        FractionalStepQuadraticIbm.h
        FractionalStepMultiphase.h
        FractionalStepMultiphase.h
        FractionalStepMultiphaseQuadraticIbm.h
        FractionalStepEulerLagrange.h
        FractionalStepMultiphaseEulerLagrange.h
        FractionalStepDirectForcing.h
        FractionalStepBoussinesq.h
        SolverFactory.h)

set(SOURCES Solver.cpp
        Poisson.cpp
        FractionalStep.cpp
        FractionalStepAxisymmetric.cpp
        FractionalStepQuadraticIbm.cpp
        FractionalStepMultiphase.cpp
        FractionalStepMultiphaseQuadraticIbm.cpp
        FractionalStepEulerLagrange.cpp
        FractionalStepMultiphaseEulerLagrange.cpp
        FractionalStepDirectForcing.cpp
        FractionalStepBoussinesq.cpp
        SolverFactory.cpp)

add_library(Solvers ${HEADERS} ${SOURCES})
target_link_libraries(Solvers Math FiniteVolumeGrid2D FiniteVolume)
install(TARGETS
        Solvers
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

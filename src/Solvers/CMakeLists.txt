set(HEADERS Solver.h
        Poisson.h
        Piso.h
        FractionalStep.h
        FractionalStepAxisymmetric.h
        FractionalStepQuadraticIbm.h
        FractionalStepMultiphase.h
        FractionalStepIncremental.h
        Coupled.h
        CoupledEquation.h
        PisoMultiphase.h
        FractionalStepIncrementalMultiphase.h
        FractionalStepMultiphase.h
        FractionalStepMultiphaseQuadraticIbm.h)

set(SOURCES Solver.cpp
        Poisson.cpp
        Piso.cpp
        FractionalStep.cpp
        FractionalStepAxisymmetric.cpp
        FractionalStepQuadraticIbm.cpp
        FractionalStepIncremental.cpp
        Coupled.cpp
        CoupledEquation.cpp
        PisoMultiphase.cpp
        FractionalStepIncrementalMultiphase.cpp
        FractionalStepMultiphase.cpp
        FractionalStepMultiphaseQuadraticIbm.cpp)

add_library(Solvers ${HEADERS} ${SOURCES})
target_link_libraries(Solvers Math FiniteVolumeGrid2D FiniteVolume)

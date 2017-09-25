set(HEADERS Solver.h
        Poisson.h
        Piso.h
        FractionalStep.h
        FractionalStepQuadraticIbm.h
        FractionalStepMultiphase.h
        FractionalStepIncremental.h
        Coupled.h
        CoupledEquation.h
        PisoMultiphase.h
        FractionalStepIncrementalMultiphase.h
        SurfaceTensionForce.h
        Celeste.h)

set(SOURCES Solver.cpp
        Poisson.cpp
        Piso.cpp
        FractionalStep.cpp
        FractionalStepQuadraticIbm.cpp
        FractionalStepIncremental.cpp
        Coupled.cpp
        CoupledEquation.cpp
        PisoMultiphase.cpp
        FractionalStepIncrementalMultiphase.cpp
        FractionalStepMultiphase.cpp
        SurfaceTensionForce.cpp
        Celeste.cpp)

add_library(Solvers ${HEADERS} ${SOURCES})

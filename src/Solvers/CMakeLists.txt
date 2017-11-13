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
        FractionalStepMultiphaseQuadraticIbm.h
        SurfaceTensionForce.h
        Celeste.h)

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
        FractionalStepMultiphaseQuadraticIbm.cpp
        SurfaceTensionForce.cpp
        Celeste.cpp
        CelesteStencil.cpp)

add_library(Solvers ${HEADERS} ${SOURCES})

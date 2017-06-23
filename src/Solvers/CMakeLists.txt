set(HEADERS Solver.h
            Poisson.h
            Piso.h
            FractionalStepSimple.h
            FractionalStepSimpleMultiphase.h
            FractionalStep.h
            Coupled.h
            CoupledEquation.h
            PisoMultiphase.h
            FractionalStepMultiphase.h
            SurfaceTensionForce.h
            ContinuumSurfaceForce.h
            Celeste.h)

set(SOURCES Solver.cpp
            Poisson.cpp
            Piso.cpp
            FractionalStepSimple.cpp
            FractionalStep.cpp
            Coupled.cpp
            CoupledEquation.cpp
            PisoMultiphase.cpp
            FractionalStepMultiphase.cpp
            FractionalStepSimpleMultiphase.cpp
            SurfaceTensionForce.cpp
            ContinuumSurfaceForce.cpp
            Celeste.cpp)

add_library(Solvers ${HEADERS} ${SOURCES})

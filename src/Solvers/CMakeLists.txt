set(HEADERS Solver.h
            Poisson.h
            Piso.h
            FractionalStep.h
            FractionalStepExperimental.h
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
            FractionalStep.cpp
            FractionalStepExperimental.cpp
            Coupled.cpp
            CoupledEquation.cpp
            PisoMultiphase.cpp
            FractionalStepMultiphase.cpp
            SurfaceTensionForce.cpp
            ContinuumSurfaceForce.cpp
            Celeste.cpp)

add_library(Solvers ${HEADERS} ${SOURCES})
target_link_libraries(Solvers Discretization ImmersedBoundary Utilities)

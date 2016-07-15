set(HEADERS PisoMultiphase.h
            FractionalStepMultiphase.h
            Cicsam.h
            Plic.h
            SurfaceTensionForce.h
            ContinuumSurfaceForce.h
            Celeste.h)

set(SOURCES PisoMultiphase.cpp
            FractionalStepMultiphase.cpp
            Cicsam.cpp
            Plic.cpp
            SurfaceTensionForce.cpp
            ContinuumSurfaceForce.cpp
            Celeste.cpp)

add_library(Multiphase ${HEADERS} ${SOURCES})
target_link_libraries(Multiphase Solvers)

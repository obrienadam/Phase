set(HEADERS Multiphase.h
            Cicsam.h
            Plic.h
            SurfaceTensionForce.h
            ContinuumSurfaceForce.h)

set(SOURCES Multiphase.cpp
            Cicsam.cpp
            Plic.cpp
            SurfaceTensionForce.cpp
            ContinuumSurfaceForce.cpp)

add_library(Multiphase ${HEADERS} ${SOURCES})
target_link_libraries(Multiphase Solvers)

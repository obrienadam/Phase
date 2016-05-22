set(HEADERS Multiphase.h
            Cicsam.h
            Plic.h
            SurfaceTensionForce.h
            ContinuumSurfaceForce.h
            HeightFunction.h)

set(SOURCES Multiphase.cpp
            Cicsam.cpp
            Plic.cpp
            SurfaceTensionForce.cpp
            ContinuumSurfaceForce.cpp
            HeightFunction.cpp)

add_library(Multiphase ${HEADERS} ${SOURCES})
target_link_libraries(Multiphase Solvers)

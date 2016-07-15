set(HEADERS ImmersedBoundary.h
            ImmersedBoundaryObject.h
            GhostCellImmersedBoundary.h
            ImmersedBoundaryContinuumSurfaceForce.h
            ImmersedBoundaryCeleste.h)

set(SOURCES ImmersedBoundary.cpp
            ImmersedBoundaryObject.cpp
            GhostCellImmersedBoundary.cpp
            ImmersedBoundaryContinuumSurfaceForce.cpp
            ImmersedBoundaryCeleste.cpp)

add_library(ImmersedBoundary ${HEADERS} ${SOURCES})
target_link_libraries(ImmersedBoundary Math Equation)

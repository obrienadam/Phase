set(HEADERS ImmersedBoundary.h
            ImmersedBoundaryObject.h
            BilinearInterpolation.cpp
            GhostCellImmersedBoundary.h
            ImmersedBoundaryContinuumSurfaceForce.h)

set(SOURCES ImmersedBoundary.cpp
            ImmersedBoundaryObject.cpp
            BilinearInterpolation.cpp
            GhostCellImmersedBoundary.cpp
            ImmersedBoundaryContinuumSurfaceForce.cpp)

add_library(ImmersedBoundary ${HEADERS} ${SOURCES})
target_link_libraries(ImmersedBoundary Multiphase)

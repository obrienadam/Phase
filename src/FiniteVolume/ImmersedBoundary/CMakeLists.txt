set(HEADERS ImmersedBoundary.h
            ImmersedBoundaryObject.h
            TranslatingImmersedBoundaryObject.h
            OscillatingImmersedBoundaryObject.h
            GhostCellStencil.h
            CutCellStencil.h
            CutCellImmersedBoundary.h
            GhostCellImmersedBoundary.h
            MovingGhostCellImmersedBoundary.h)

set(SOURCES ImmersedBoundary.cpp
            ImmersedBoundaryObject.cpp
            TranslatingImmersedBoundaryObject.cpp
            OscillatingImmersedBoundaryObject.cpp
            GhostCellStencil.cpp
            CutCellStencil.cpp
            CutCellImmersedBoundary.cpp
            GhostCellImmersedBoundary.cpp
            MovingGhostCellImmersedBoundary.cpp)

add_library(ImmersedBoundary ${HEADERS} ${SOURCES})
target_link_libraries(ImmersedBoundary Math Equation)

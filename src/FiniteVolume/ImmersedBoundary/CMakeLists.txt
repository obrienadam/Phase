set(HEADERS ImmersedBoundary.h
            ImmersedBoundaryObject.h
            GhostCellImmersedBoundaryObject.h
            ForcingCellImmersedBoundaryObject.h
            LeastSquaresImmersedBoundaryObject.h
            GhostCellStencil.h
            CutCell.h
            CutCellLink.h
            ImmersedBoundaryStencil.h
            ForcingCellStencil.h
            LeastSquaresStencil.h
            GhostCellImmersedBoundary.h
            ForcingCellImmersedBoundary.h
            MovingGhostCellImmersedBoundary.h)

set(SOURCES ImmersedBoundary.cpp
            ImmersedBoundaryObject.cpp
            GhostCellImmersedBoundaryObject.cpp
            ForcingCellImmersedBoundaryObject.cpp
            LeastSquaresImmersedBoundaryObject.cpp
            ImmersedBoundaryStencil.cpp
            ForcingCellStencil.cpp
            LeastSquaresStencil.cpp
            GhostCellStencil.cpp
            CutCell.cpp
            CutCellLink.cpp
            GhostCellImmersedBoundary.cpp
            MovingGhostCellImmersedBoundary.cpp)

add_library(ImmersedBoundary ${HEADERS} ${SOURCES})
target_link_libraries(ImmersedBoundary Math Equation)

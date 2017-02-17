set(HEADERS ImmersedBoundary.h
            ImmersedBoundaryObject.h
            TranslatingImmersedBoundaryObject.h
            GhostCellImmersedBoundary.h)

set(SOURCES ImmersedBoundary.cpp
            ImmersedBoundaryObject.cpp
            TranslatingImmersedBoundaryObject.cpp
            GhostCellImmersedBoundary.cpp)

add_library(ImmersedBoundary ${HEADERS} ${SOURCES})
target_link_libraries(ImmersedBoundary Math Equation)

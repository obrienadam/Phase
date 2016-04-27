set(HEADERS ImmersedBoundary.h
            ImmersedBoundaryObject.h
            BilinearInterpolation.cpp)

set(SOURCES ImmersedBoundary.cpp
            ImmersedBoundaryObject.cpp
            BilinearInterpolation.cpp)

add_library(ImmersedBoundary ${HEADERS} ${SOURCES})
target_link_libraries(ImmersedBoundary Multiphase)

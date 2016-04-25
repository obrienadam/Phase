set(HEADERS ImmersedBoundary.h
            BilinearInterpolation.cpp)

set(SOURCES ImmersedBoundary.cpp
            BilinearInterpolation.cpp)

add_library(ImmersedBoundary ${HEADERS} ${SOURCES})
target_link_libraries(ImmersedBoundary Multiphase)

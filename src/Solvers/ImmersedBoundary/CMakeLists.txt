set(HEADERS ImmersedBoundary.h)

set(SOURCES ImmersedBoundary.cpp)

add_library(ImmersedBoundary ${HEADERS} ${SOURCES})
target_link_libraries(ImmersedBoundary Multiphase)

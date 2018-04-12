include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB HEADERS
        StructuredGrid2D/*.h
        FiniteDifference/*.h)

file(GLOB SOURCES
        StructuredGrid2D/*.cpp)

add_library(phase_2d_structured ${HEADERS} ${SOURCES})
target_link_libraries(phase_2d_structured phase_2d_geometry)
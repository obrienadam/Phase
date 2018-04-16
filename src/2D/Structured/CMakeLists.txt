include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB HEADERS
        StructuredGrid2D/*.h
        Solvers/*.h
        FiniteDifference/*.h)

file(GLOB SOURCES
        StructuredGrid2D/*.cpp
        Solvers/*.cpp)

add_library(phase_2d_structured ${HEADERS} ${SOURCES} utilities/Phase2DStructuredPartitionGrid.cpp Solvers/StructuredFractionalStep.cpp Solvers/StructuredFractionalStep.h)
target_link_libraries(phase_2d_structured phase_2d_geometry)

add_executable(phase-2d-structured-partition-grid utilities/Phase2DStructuredPartitionGrid.cpp)
target_link_libraries(phase-2d-structured-partition-grid phase_system)
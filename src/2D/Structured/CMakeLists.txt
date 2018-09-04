include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB HEADERS
        StructuredGrid2D/*.h
        Solvers/*.h
        FiniteVolume/*/*.h
        PostProcessing/*.h)

file(GLOB SOURCES
        StructuredGrid2D/*.cpp
        StructuredGrid2D/*.tpp
        Solvers/*.cpp
        FiniteVolume/*/*.cpp
        FiniteVolume/*/*.tpp
        PostProcessing/*.cpp)

add_library(phase_2d_structured ${HEADERS} ${SOURCES})
target_link_libraries(phase_2d_structured phase_2d_geometry phase_math cgns metis)

add_executable(phase-2d-structured-partition-grid utilities/Phase2DStructuredPartitionGrid.cpp)
target_link_libraries(phase-2d-structured-partition-grid phase_system)

add_executable(phase-2d-structured modules/Phase2DStructured)
target_link_libraries(phase-2d-structured phase_system phase_2d_structured)

install(TARGETS
        phase_2d_structured
        phase-2d-structured-partition-grid
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

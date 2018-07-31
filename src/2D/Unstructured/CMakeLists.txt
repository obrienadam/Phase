include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB HEADERS
        FiniteVolume/*/*.h
        FiniteVolumeGrid2D/*/*.h
        FiniteVolumeGrid2D/*.h
        Solvers/*.h
        PostProcessing/*.h)

file(GLOB SOURCES
        FiniteVolume/*/*.cpp
        FiniteVolume/*/*.tpp
        FiniteVolumeGrid2D/*/*.cpp
        FiniteVolumeGrid2D/*.cpp
        FiniteVolumeGrid2D/*.tpp
        Solvers/*.cpp
        PostProcessing/*.cpp)

add_library(phase_2d_unstructured ${HEADERS} ${SOURCES})
target_link_libraries(phase_2d_unstructured phase_2d_geometry phase_math cgns metis)

add_executable(phase-2d-unstructured modules/Phase2DUnstructured.cpp)
target_link_libraries(phase-2d-unstructured phase_2d_unstructured)

add_executable(phase-2d-unstructured-partition-grid utilities/PhasePartitionGrid.cpp)
target_link_libraries(phase-2d-unstructured-partition-grid phase_system phase_2d_unstructured)

install(TARGETS
        phase_2d_unstructured
        phase-2d-unstructured
        phase-2d-unstructured-partition-grid
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

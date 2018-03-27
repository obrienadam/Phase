include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB HEADERS
        FiniteVolume/*/*.h
        FiniteVolumeGrid2D/*.h
        FiniteVolumeGrid2D/*/*.h
        Solvers/*.h
        PostProcessing/*.h)

file(GLOB SOURCES
        FiniteVolume/*/*.cpp
        FiniteVolumeGrid2D/*.cpp
        FiniteVolumeGrid2D/*/*.cpp
        Solvers/*.cpp
        PostProcessing/*.cpp)

add_library(phase_2d_unstructured ${HEADERS} ${SOURCES} PostProcessing/PostProcessing.cpp PostProcessing/PostProcessing.h)
target_link_libraries(phase_2d_unstructured phase_2d_geometry phase_math cgns metis)

add_executable(phase-2d-unstructured modules/Phase2DUnstructured.cpp ../Structured/StructuredGrid2D/StructuredGrid2D.cpp ../Structured/StructuredGrid2D/StructuredGrid2D.h)
target_link_libraries(phase-2d-unstructured phase_2d_unstructured)

install(TARGETS
        phase_2d_unstructured
        phase-2d-unstructured
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
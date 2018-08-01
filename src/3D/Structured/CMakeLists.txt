file(GLOB HEADERS */*.h */*/*.h)
file(GLOB SOURCES */*.cpp */*.tpp */*/*.cpp */*/*.tpp)

add_library(phase_3d_structured ${HEADERS} ${SOURCES})
target_link_libraries(phase_3d_structured phase_3d_geometry phase_math)

add_executable(phase-3d-structured modules/Phase3DStructured.cpp)
target_link_libraries(phase-3d-structured phase_3d_structured)

install(TARGETS
    phase_3d_structured
    phase-3d-structured
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)

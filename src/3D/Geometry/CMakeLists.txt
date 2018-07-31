file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)

add_library(phase_3d_geometry ${HEADERS} ${SOURCES})

install(TARGETS
        phase_3d_geometry
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

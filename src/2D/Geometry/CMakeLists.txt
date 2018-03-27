file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)

add_library(phase_2d_geometry ${HEADERS} ${SOURCES})

install(TARGETS
        phase_2d_geometry
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

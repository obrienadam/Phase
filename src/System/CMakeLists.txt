set(HEADERS Input.h
            CommandLine.h
            Exception.h)

set(SOURCES Input.cpp
            CommandLine.cpp
            Exception.cpp)

add_library(System ${HEADERS} ${SOURCES})
target_link_libraries(System boost_program_options boost_system boost_filesystem)

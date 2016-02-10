set(HEADERS Input.h
            CommandLine.h
            Exception.h)
set(SOURCES Input.cpp
            CommandLine.cpp
            Exception.cpp)

add_library(InputOutput ${HEADERS} ${SOURCES})
target_link_libraries(InputOutput boost_program_options)

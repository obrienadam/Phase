set(HEADERS Input.h
            CommandLine.h
            Exception.h
            RunControl.h)

set(SOURCES Input.cpp
            CommandLine.cpp
            Exception.cpp
            RunControl.cpp)

add_library(System ${HEADERS} ${SOURCES})
target_link_libraries(System)

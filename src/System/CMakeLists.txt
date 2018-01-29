set(HEADERS Input.h
            CommandLine.h
            Exception.h
            Time.h
            RunControl.h)

set(SOURCES Input.cpp
            CommandLine.cpp
            Exception.cpp
            Time.cpp
            RunControl.cpp)

add_library(System ${HEADERS} ${SOURCES})
target_link_libraries(System PostProcessing ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
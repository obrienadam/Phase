set(HEADERS CrankNicolson.h
            AdamsBashforth.h)

set(SOURCES CrankNicolson.cpp
            AdamsBashforth.cpp)

add_library(Discretization ${HEADERS} ${SOURCES})
target_link_libraries(Discretization Field Math)

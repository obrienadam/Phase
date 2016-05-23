set(HEADERS Equation.h)

set(SOURCES Equation.cpp
            Equation.tpp)

add_library(Equation ${HEADERS} ${SOURCES})
target_link_libraries(Equation Field Math)

set(HEADERS Equation.h
            FiniteVolumeEquation.h
            TimeDerivative.h
            Divergence.h
            Laplacian.h
            Source.h)

set(SOURCES Equation.cpp
            Equation.tpp
            TimeDerivative.cpp
            Divergence.cpp
            Laplacian.cpp
            Source.cpp)

add_library(Equation ${HEADERS} ${SOURCES})
target_link_libraries(Equation Field Math)

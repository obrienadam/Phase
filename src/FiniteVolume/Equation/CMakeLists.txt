set(HEADERS Equation.h
            Matrix.h
            SparseMatrix.h
            SparseVector.h
            Diffusion.h)

set(SOURCES Equation.cpp
            Matrix.cpp
            SparseMatrix.cpp
            SparseVector.cpp
            Diffusion.cpp)

add_library(Equation ${HEADERS} ${SOURCES})

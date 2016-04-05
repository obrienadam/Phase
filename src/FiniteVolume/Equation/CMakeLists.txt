set(HEADERS Equation.h
            Matrix.h
            SparseMatrix.h
            SparseVector.h)

set(SOURCES Equation.cpp
            Equation.tpp
            Matrix.cpp
            SparseMatrix.cpp
            SparseVector.cpp)

add_library(Equation ${HEADERS} ${SOURCES})
target_link_libraries(Equation Field)

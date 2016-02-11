set(HEADERS Equation.h
            Matrix.h
            SparseMatrix.h
            SparseVector.h
            Term.h
            AdvectionTerm.h
            DiffusionTerm.h)

set(SOURCES Equation.cpp
            Matrix.cpp
            SparseMatrix.cpp
            SparseVector.cpp
            Term.cpp
            AdvectionTerm.cpp
            DiffusionTerm.cpp)

add_library(Equation ${HEADERS} ${SOURCES})

set(HEADERS Equation.h
            Matrix.h
            SparseMatrix.h
            SparseVector.h
            Term.h
            AdvectionTerm.h
            DiffusionTerm.h
            TemporalTerm.h)

set(SOURCES Equation.cpp
            Equation.tpp
            Matrix.cpp
            SparseMatrix.cpp
            SparseVector.cpp
            Term.cpp
            AdvectionTerm.cpp
            DiffusionTerm.cpp
            TemporalTerm.cpp)

add_library(Equation ${HEADERS} ${SOURCES})
target_link_libraries(Equation Field)

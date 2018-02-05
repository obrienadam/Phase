set(HEADERS StaticMatrix.h
        Matrix.h
        StaticMatrix.h
        SparseMatrixSolver.h
        EigenSparseMatrixSolver.h
        TrilinosBelosSparseMatrixSolver.h
        TrilinosMueluSparseMatrixSolver.h
        Vector.h
        Algorithm.h)

set(SOURCES Matrix.cpp
        SparseMatrixSolver.cpp
        EigenSparseMatrixSolver.cpp
        TrilinosBelosSparseMatrixSolver.cpp
        TrilinosMueluSparseMatrixSolver.cpp
        Vector.cpp)

add_library(Math ${HEADERS} ${SOURCES})
target_link_libraries(Math Communicator
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        lapacke
        ${Trilinos_LIBRARIES}
        ${Trilinos_TPL_LIBRARIES})

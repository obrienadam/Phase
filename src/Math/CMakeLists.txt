set(HEADERS Matrix.h
            SparseMatrixSolver.h
            EigenSparseMatrixSolver.h
            Vector.h
            HypreSparseMatrixSolver.h
            PetscSparseMatrixSolver.h
            BisectionSearch.h
            SecantSearch.h
            BilinearInterpolation.h)

set(SOURCES Matrix.cpp
            SparseMatrixSolver.cpp
            EigenSparseMatrixSolver.cpp
            Vector.cpp
            HypreSparseMatrixSolver.cpp
            PetscSparseMatrixSolver.cpp
            BisectionSearch.cpp
            SecantSearch.cpp
            BilinearInterpolation.cpp)

add_library(Math ${HEADERS} ${SOURCES})
target_link_libraries(Math lapacke lapack openblas HYPRE petsc)

set(HEADERS StaticMatrix.h
        Matrix.h
        BlockMatrix.h
        StaticMatrix.h
        SparseMatrixSolver.h
        EigenSparseMatrixSolver.h
        TrilinosBelosSparseMatrixSolver.h
        #TrilinosMueluSparseMatrixSolver.h
        Vector.h
        Algorithm.h
        Interpolation.h
        LinearInterpolation.h
        BilinearInterpolation.h
        QuadraticNormalInterpolation.h
        InterpolationFunction.h)

set(SOURCES Matrix.cpp
        BlockMatrix.cpp
        SparseMatrixSolver.cpp
        EigenSparseMatrixSolver.cpp
        TrilinosBelosSparseMatrixSolver.cpp
        #TrilinosMueluSparseMatrixSolver.cpp
        Vector.cpp
        LinearInterpolation.cpp
        BilinearInterpolation.cpp
        QuadraticNormalInterpolation.cpp
        InterpolationFunction.cpp Algorithm.h)

add_library(Math ${HEADERS} ${SOURCES})

target_link_libraries(Math ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        ${Trilinos_LIBRARIES}
        ${Trilinos_TPL_LIBRARIES})

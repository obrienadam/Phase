set(HEADERS Matrix.h
    SparseMatrixSolver.h
    EigenSparseMatrixSolver.h
    TrilinosSparseMatrixSolver.h
    Vector.h
    BisectionSearch.h
    SecantSearch.h
    Interpolation.h
    LinearInterpolation.h
    BilinearInterpolation.h
    QuadraticNormalInterpolation.h
    InterpolationFunction.h)

set(SOURCES Matrix.cpp
    SparseMatrixSolver.cpp
    EigenSparseMatrixSolver.cpp
    TrilinosSparseMatrixSolver.cpp
    Vector.cpp
    BisectionSearch.cpp
    SecantSearch.cpp
    LinearInterpolation.cpp
    BilinearInterpolation.cpp
    QuadraticNormalInterpolation.cpp
    InterpolationFunction.cpp)

add_library(Math ${HEADERS} ${SOURCES})

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    target_link_libraries(Math ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

target_link_libraries(Math ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

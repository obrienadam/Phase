set(HEADERS Matrix.h
    SparseMatrixSolver.h
    Vector.h
    BisectionSearch.h
    SecantSearch.h
    Interpolation.h
    BilinearInterpolation.h
    QuadraticNormalInterpolation.h)

set(SOURCES Matrix.cpp
    SparseMatrixSolver.cpp
    Vector.cpp
    BisectionSearch.cpp
    SecantSearch.cpp
    BilinearInterpolation.cpp
    QuadraticNormalInterpolation.cpp)

if(PHASE_USE_EIGEN)
    set(HEADERS ${HEADERS} EigenSparseMatrixSolver.h)
    set(SOURCES ${SOURCES} EigenSparseMatrixSolver.cpp)
endif()

if(PHASE_USE_HYPRE)
    set(HEADERS ${HEADERS} HypreSparseMatrixSolver.h)
    set(SOURCES ${SOURCES} HypreSparseMatrixSolver.cpp)
endif()

if(PHASE_USE_PETSC)
    set(HEADERS ${HEADERS} PetscSparseMatrixSolver.h)
    set(SOURCES ${SOURCES} PetscSparseMatrixSolver.cpp)
endif()

if(PHASE_USE_TRILINOS)
    set(HEADERS ${HEADERS} TrilinosSparseMatrixSolver.h)
    set(SOURCES ${SOURCES} TrilinosSparseMatrixSolver.cpp)
endif()

add_library(Math ${HEADERS} ${SOURCES})

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    target_link_libraries(Math lapacke lapack openblas)
endif()

if(PHASE_USE_PETSC)
    target_link_libraries(Math petsc)
endif()

if(PHASE_USE_HYPRE)
    target_link_libraries(Math HYPRE)
endif()

if(PHASE_USE_TRILINOS)
    target_link_libraries(Math
        teuchoscore
        teuchosparameterlist
        teuchoscomm
        teuchoskokkoscompat
        teuchosnumerics
        teuchosremainder
        kokkoscore
        tpetra
        tpetrakernels
        belos
        ifpack2
        amesos2
        zoltan2)
endif()

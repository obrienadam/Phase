set(HEADERS StaticMatrix.h
        Matrix.h
        StaticMatrix.h
        SparseMatrixSolver.h
        EigenSparseMatrixSolver.h
        TrilinosSparseMatrixSolver.h
        TrilinosBelosSparseMatrixSolver.h
        TrilinosAmesosSparseMatrixSolver.h
        TrilinosMueluSparseMatrixSolver.h
        SparseMatrixSolverFactory.h
        Equation.h
        Vector.h
        Algorithm.h)

set(SOURCES Matrix.cpp
        SparseMatrixSolver.cpp
        EigenSparseMatrixSolver.cpp
        TrilinosSparseMatrixSolver.cpp
        TrilinosBelosSparseMatrixSolver.cpp
        TrilinosAmesosSparseMatrixSolver.cpp
        TrilinosMueluSparseMatrixSolver.cpp
        SparseMatrixSolverFactory.cpp
        Equation.cpp
        Vector.cpp)

add_library(phase_math ${HEADERS} ${SOURCES})

target_link_libraries(phase_math
        phase_system
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        #lapacke
        ${Trilinos_LIBRARIES}
        ${Trilinos_TPL_LIBRARIES})

install(TARGETS
        phase_math
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

#ifndef HYPRE_SPARSE_MATRIX_SOLVER_H
#define HYPRE_SPARSE_MATRIX_SOLVER_H

#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>

#include "Communicator.h"
#include "Equation.h"

class HypreSparseMatrixSolver
{
public:
    HypreSparseMatrixSolver(const Communicator &comm);

    void set(int i, int j, Scalar value);

private:
    const Communicator& comm_;
    HYPRE_IJMatrix ijMatrix_;
    HYPRE_ParCSRMatrix csrMatrix_;
};

#endif

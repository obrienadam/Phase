#include "SparseMatrixSolver.h"

Scalar SparseMatrixSolver::solve(const Vector &x0)
{
    return solve();
}

void SparseMatrixSolver::setMaxPreconditionerUses(int maxPreconditionerUses)
{
    maxPreconUses_ = maxPreconditionerUses;
    nPreconUses_ = maxPreconUses_; // done to ensure the preconditioner will be recomputed
}

void SparseMatrixSolver::printStatus(const std::string &msg) const
{
    printf("%s iterations = %d, error = %lf.\n", msg.c_str(), nIters(), error());
}

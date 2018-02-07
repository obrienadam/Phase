#include "SparseMatrixSolver.h"

Scalar SparseMatrixSolver::solve(const Vector &x0)
{
    setGuess(x0);
    return solve();
}

void SparseMatrixSolver::printStatus(const std::string &msg) const
{
    printf("%s iterations = %d, error = %lf.\n", msg.c_str(), nIters(), error());
}

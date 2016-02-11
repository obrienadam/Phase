#include "SparseMatrix.h"

SparseMatrix::SparseMatrix(int nRows, int nCols, int nnz)
    :
      Eigen::SparseMatrix<Scalar>(nRows, nCols)
{
    reserve(nRows*nnz);
}

void SparseMatrix::assemble(const std::vector<Eigen::Triplet<Scalar> > &entries)
{
    setFromTriplets(entries.begin(), entries.end());
    solver_.compute(*this);
}

SparseVector SparseMatrix::solve(const SparseVector &b) const
{
    return solver_.solve(b);
}

SparseVector SparseMatrix::solve(const SparseVector &b, const SparseVector &x0) const
{
    return solver_.solveWithGuess(b, x0);
}

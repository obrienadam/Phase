#include "SparseMatrix.h"
#include <stdio.h>

SparseMatrix::SparseMatrix(int nRows, int nCols, int nnz)
    :
      Eigen::SparseMatrix<Scalar>(nRows, nCols)
{
    reserve(nRows*nnz);
    setZero();
    solver_.setMaxIterations(25);
    solver_.setTolerance(1e-12);
}

SparseMatrix::SparseMatrix(const SparseMatrix &other)
    :
      Eigen::SparseMatrix<Scalar>(other)
{
    solver_.setMaxIterations(other.solver_.maxIterations());
    solver_.setTolerance(other.solver_.tolerance());
}

SparseMatrix& SparseMatrix::operator =(const SparseMatrix& rhs)
{
    Eigen::SparseMatrix<Scalar>::operator =(rhs);

    return *this;
}

void SparseMatrix::assemble(const std::vector<Eigen::Triplet<Scalar> > &entries)
{
    setFromTriplets(entries.begin(), entries.end());
}

SparseVector SparseMatrix::solve(const SparseVector &b) const
{
    solver_.compute(*this);
    return solver_.solve(b);
}

SparseVector SparseMatrix::solve(const SparseVector &b, const SparseVector &x0) const
{
    solver_.compute(*this);
    return solver_.solveWithGuess(b, x0);
}

#include "SparseMatrix.h"

SparseMatrix::SparseMatrix(int nRows, int nCols, int nnz)
    :
      Eigen::SparseMatrix<Scalar>(nRows, nCols)
{
    reserve(nRows*nnz);
}

SparseMatrix::SparseMatrix(const SparseMatrix &other)
    :
      Eigen::SparseMatrix<Scalar>(other)
{

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
    return solver_.solveWithGuess(b, x0);
}

#include "SparseMatrix.h"

SparseMatrix::SparseMatrix(int nRows, int nCols, int nnz)
    :
      Eigen::SparseMatrix<Scalar>(nRows, nCols)
{
    reserve(nRows*nnz);
    entries_.reserve(nRows*nnz);
}

void SparseMatrix::addEntry(int row, int col, Scalar entry)
{
    entries_.push_back(Eigen::Triplet<Scalar>(row, col, entry));
}

void SparseMatrix::addEntries(const std::vector<int> &rows, const std::vector<int> &cols, const std::vector<Scalar> &entries)
{
    for(int i = 0, end = entries.size(); i < end; ++i)
        entries_.push_back(Eigen::Triplet<Scalar>(rows[i], cols[i], entries[i]));
}

void SparseMatrix::assemble()
{
    setFromTriplets(entries_.begin(), entries_.end());
    solver_.compute(*this);
    entries_.clear();
}

SparseVector SparseMatrix::solve(const SparseVector &b) const
{
    return solver_.solve(b);
}

SparseVector SparseMatrix::solve(const SparseVector &b, const SparseVector &x0) const
{
    return solver_.solveWithGuess(b, x0);
}

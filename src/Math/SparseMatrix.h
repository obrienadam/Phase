#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <vector>

#include <eigen3/Eigen/Sparse>

#include "Types.h"
#include "SparseVector.h"

class SparseMatrix : public Eigen::SparseMatrix<Scalar>
{
public:

    typedef Eigen::Triplet<Scalar> Triplet;

    enum Preconditioner{IncompleteLUT, NoPreconditioner};

    SparseMatrix() {}
    SparseMatrix(int nRows, int nCols, int nnz);
    SparseMatrix(const SparseMatrix& other);

    SparseMatrix& operator=(const SparseMatrix& rhs);

    void addEntry(int row, int col, Scalar entry);
    void addEntries(const std::vector<int>& rows, const std::vector<int>& cols, const std::vector<Scalar>& entries);

    void assemble(const std::vector< Eigen::Triplet<Scalar> >& entries);

    SparseVector solve(const SparseVector& b, Preconditioner precon = IncompleteLUT) const;//SparseVector solve(const SparseVector &b, const SparseVector &x0, Preconditioner precon = IncompleteLUT) const;

    Scalar error() const { return error_; }
    int nIterations() const { return nIters_; }

private:

    mutable Scalar error_;
    mutable size_t nIters_;

    mutable Eigen::BiCGSTAB< SparseMatrix, Eigen::IncompleteLUT<Scalar> > solverIncompleteLUT_;
    mutable Eigen::BiCGSTAB< SparseMatrix, Eigen::IdentityPreconditioner > solverNoPreconditioner_;
};

#endif

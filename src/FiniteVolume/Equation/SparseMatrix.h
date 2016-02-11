#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <vector>

#include <eigen3/Eigen/Sparse>

#include "Types.h"
#include "SparseVector.h"

class SparseMatrix : public Eigen::SparseMatrix<Scalar>
{
public:

    SparseMatrix(int nRows, int nCols, int nnz);

    void addEntry(int row, int col, Scalar entry);
    void addEntries(const std::vector<int>& rows, const std::vector<int>& cols, const std::vector<Scalar>& entries);

    void assemble(const std::vector< Eigen::Triplet<Scalar> >& entries);

    SparseVector solve(const SparseVector& b) const;
    SparseVector solve(const SparseVector &b, const SparseVector &x0) const;

    Scalar error() const { return solver_.error(); }
    int nIterations() const { return solver_.iterations(); }

private:

    Eigen::BiCGSTAB< SparseMatrix, Eigen::IncompleteLUT<Scalar> > solver_;
};

#endif

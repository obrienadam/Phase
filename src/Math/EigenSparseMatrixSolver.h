#ifndef PHASE_EIGEN_SPARSE_MATRIX_SOLVER_H
#define PHASE_EIGEN_SPARSE_MATRIX_SOLVER_H

#include <vector>

#include <eigen3/Eigen/SparseLU>

#include "SparseMatrixSolver.h"

class EigenSparseMatrixSolver : public SparseMatrixSolver
{
public:

    typedef Eigen::Triplet<Scalar> Triplet;
    typedef Eigen::SparseMatrix<Scalar> EigenSparseMatrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> EigenVector;
    typedef Eigen::SparseLU<EigenSparseMatrix> SparseLUSolver;

    EigenSparseMatrixSolver();

    Type type() const
    { return EIGEN; }

    void setRank(int rank);

    void setRank(int rowRank, int colRank);

    void set(const CoefficientList &coeffs);

    void setGuess(const Vector &x0);

    void setRhs(const Vector &rhs);

    Scalar solve();

    Scalar solve(const Vector &x0);

    Scalar x(Index idx) const
    { return x_[idx]; }

    int nIters() const
    { return 1; }

    Scalar error() const
    { return 0.; }

    bool supportsMPI() const
    { return false; }

    std::shared_ptr<SparseMatrixSolver> newSparseMatrixSolver() const
    { return std::shared_ptr<EigenSparseMatrixSolver>(new EigenSparseMatrixSolver()); }

private:

    EigenSparseMatrix mat_;

    EigenVector x_, rhs_;

    SparseLUSolver solver_;
};

#endif

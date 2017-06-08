#ifndef EIGEN_SPARSE_MATRIX_SOLVER_H
#define EIGEN_SPARSE_MATRIX_SOLVER_H

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

    void setRank(int rank);

    void set(const CoefficientList &coeffs);

    void setGuess(const Vector &x0);

    void setRhs(const Vector &rhs);

    Scalar solve();

    Scalar solve(const Vector &x0);

    void mapSolution(ScalarFiniteVolumeField &field);

    void mapSolution(VectorFiniteVolumeField &field);

    void setMaxIters(int maxIters)
    { }

    void setToler(Scalar toler)
    { }

    void setDropToler(Scalar toler)
    { }

    void setFillFactor(int fill)
    { }

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

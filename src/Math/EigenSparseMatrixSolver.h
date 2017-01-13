#ifndef EIGEN_SPARSE_MATRIX_SOLVER_H
#define EIGEN_SPARSE_MATRIX_SOLVER_H

#include <vector>

#include <eigen3/Eigen/Sparse>

#include "SparseMatrixSolver.h"

class EigenSparseMatrixSolver : public SparseMatrixSolver
{
public:

    typedef Eigen::Triplet<Scalar> Triplet;
    typedef Eigen::SparseMatrix<Scalar> EigenSparseMatrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> EigenVector;
    typedef Eigen::Map<EigenVector> EigenVectorMap;
    typedef Eigen::BiCGSTAB< EigenSparseMatrix, Eigen::IncompleteLUT<Scalar> > BiCGSTABIncompleteLUT;
    typedef Eigen::BiCGSTAB< EigenSparseMatrix, Eigen::IdentityPreconditioner > BiCGSTABIdentityPreconditioner;
    typedef Eigen::BiCGSTAB< EigenSparseMatrix, Eigen::DiagonalPreconditioner<Scalar> > BiCGSTABDiagonalPreconditioner;

    EigenSparseMatrixSolver();

    void setRank(int rank);
    void set(const CoefficientList& coeffs);
    void setRhs(const Vector& rhs);

    Scalar solve();

    virtual void mapSolution(ScalarFiniteVolumeField& field);
    virtual void mapSolution(VectorFiniteVolumeField& field);

    virtual void setMaxIters(int maxIters) { solver_.setMaxIterations(maxIters); }
    virtual void setToler(Scalar toler) { solver_.setTolerance(toler); }
    virtual void setFillFactor(int fill) { solver_.preconditioner().setFillfactor(fill); }

    virtual int nIters() const { return solver_.iterations(); }
    virtual Scalar error() const { return solver_.error(); }

    std::shared_ptr<SparseMatrixSolver> newSparseMatrixSolver() const { return std::make_shared<EigenSparseMatrixSolver>(); }

private:

    EigenSparseMatrix mat_;
    EigenVector x_, rhs_;
    BiCGSTABIncompleteLUT solver_;
};

#endif

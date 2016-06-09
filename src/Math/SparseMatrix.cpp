#include "SparseMatrix.h"
#include "Exception.h"

SparseMatrix::SparseMatrix(int nRows, int nCols, int nnz)
    :
      Eigen::SparseMatrix<Scalar>(nRows, nCols)
{
    reserve(nRows*nnz);
    solverNoPreconditioner_.setTolerance(1e-12);
    solverNoPreconditioner_.setMaxIterations(50);
    solverIncompleteLUT_.setTolerance(1e-12);
    solverIncompleteLUT_.setMaxIterations(50);
}

SparseMatrix::SparseMatrix(const SparseMatrix &other)
    :
      Eigen::SparseMatrix<Scalar>(other)
{
    solverNoPreconditioner_.setTolerance(other.solverNoPreconditioner_.tolerance());
    solverNoPreconditioner_.setMaxIterations(other.solverNoPreconditioner_.maxIterations());
    solverIncompleteLUT_.setTolerance(other.solverIncompleteLUT_.tolerance());
    solverIncompleteLUT_.setMaxIterations(other.solverIncompleteLUT_.maxIterations());
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

SparseVector SparseMatrix::solve(const SparseVector &b, Preconditioner precon) const
{
    SparseVector x;
    Eigen::ComputationInfo info;

    switch(precon)
    {
    case NoPreconditioner:
        solverNoPreconditioner_.compute(*this);
        x = solverNoPreconditioner_.solve(b);
        error_ = solverNoPreconditioner_.error();
        nIters_ = solverNoPreconditioner_.iterations();
        info = solverNoPreconditioner_.info();
        break;

    case IncompleteLUT:
        solverIncompleteLUT_.compute(*this);
        x = solverIncompleteLUT_.solve(b);
        error_ = solverIncompleteLUT_.error();
        nIters_ = solverIncompleteLUT_.iterations();
        info = solverIncompleteLUT_.info();
        break;
    }

    //if(info == Eigen::NoConvergence)
    //    throw Exception("SparseMatrix", "solve", "iterative procedure did not converge.");

    return x;
}

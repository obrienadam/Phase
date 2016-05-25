#include "SparseMatrix.h"
#include <stdio.h>

SparseMatrix::SparseMatrix(int nRows, int nCols, int nnz)
    :
      Eigen::SparseMatrix<Scalar>(nRows, nCols)
{
    reserve(nRows*nnz);
    setZero();
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
    switch(precon)
    {
    case NoPreconditioner:
        solverNoPreconditioner_.compute(*this);
        error_ = solverNoPreconditioner_.error();
        nIters_ = solverNoPreconditioner_.iterations();
        return solverNoPreconditioner_.solve(b);

    case IncompleteLUT:
        solverIncompleteLUT_.compute(*this);
        error_ = solverIncompleteLUT_.error();
        nIters_ = solverIncompleteLUT_.iterations();
        return solverIncompleteLUT_.solve(b);

    }
}

SparseVector SparseMatrix::solve(const SparseVector &b, const SparseVector &x0, Preconditioner precon) const
{
    switch(precon)
    {
    case NoPreconditioner:
        solverNoPreconditioner_.compute(*this);
        error_ = solverNoPreconditioner_.error();
        nIters_ = solverNoPreconditioner_.iterations();
        return solverNoPreconditioner_.solveWithGuess(b, x0);

    case IncompleteLUT:
        solverIncompleteLUT_.compute(*this);
        error_ = solverIncompleteLUT_.error();
        nIters_ = solverIncompleteLUT_.iterations();
        return solverIncompleteLUT_.solveWithGuess(b, x0);

    }
}

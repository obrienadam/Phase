#include "EigenSparseMatrixSolver.h"

EigenSparseMatrixSolver::EigenSparseMatrixSolver()
{

}

void EigenSparseMatrixSolver::setRank(int rank)
{
    mat_.resize(rank, rank);
    x_.resize(rank);
    rhs_.resize(rank);
}

void EigenSparseMatrixSolver::set(const CoefficientList &coeffs)
{
    std::vector<Triplet> triplets;
    triplets.reserve(5*coeffs.size());

    for(int i = 0, end = coeffs.size(); i < end; ++i)
        for(const auto& entry: coeffs[i])
            triplets.push_back(Triplet(i, entry.first, entry.second));

    mat_.setFromTriplets(triplets.begin(), triplets.end());
}

void EigenSparseMatrixSolver::setGuess(const Vector &x0)
{
    for(int i = 0, end = x0.size(); i < end; ++i)
        x_(i) = x0(i);
}

void EigenSparseMatrixSolver::setRhs(const Vector &rhs)
{
    for(int i = 0, end = rhs.size(); i < end; ++i)
        rhs_(i) = rhs(i);
}

Scalar EigenSparseMatrixSolver::solve()
{
    if(nPreconUses_ == maxPreconUses_)
    {
        solver_.compute(mat_);
        nPreconUses_ = 0;
    }
    else
        ++nPreconUses_;

    x_ = solver_.solve(rhs_);
    return solver_.error();
}

Scalar EigenSparseMatrixSolver::solve(const Vector &x0)
{
    if(nPreconUses_ >= maxPreconUses_)
    {
        solver_.compute(mat_);
        nPreconUses_ = 1;
    }
    else
        ++nPreconUses_;

    for(int i = 0, end = x0.size(); i < end; ++i)
        x_(i) = x0(i);

    x_ = solver_.solveWithGuess(rhs_, x_);
    return solver_.error();
}

void EigenSparseMatrixSolver::mapSolution(ScalarFiniteVolumeField &field)
{
    for(const Cell& cell: field.grid.localActiveCells())
        field(cell) = x_[cell.index(0)];
}

void EigenSparseMatrixSolver::mapSolution(VectorFiniteVolumeField &field)
{
    const Size nActiveCells = field.grid.nLocalActiveCells();
    for(const Cell& cell: field.grid.localActiveCells())
    {
        Vector2D& vec = field(cell);
        vec.x = x_[cell.index(0)];
        vec.y = x_[cell.index(0) + nActiveCells];
    }
}

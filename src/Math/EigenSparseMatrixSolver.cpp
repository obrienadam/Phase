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

void EigenSparseMatrixSolver::setRhs(const Vector &rhs)
{
    for(int i = 0, end = rhs.size(); i < end; ++i)
        rhs_(i) = rhs(i);
}

Scalar EigenSparseMatrixSolver::solve()
{
    solver_.compute(mat_);
    x_ = solver_.solve(rhs_);
    return solver_.error();
}

void EigenSparseMatrixSolver::mapSolution(ScalarFiniteVolumeField &field)
{
    for(const Cell& cell: field.grid.activeCells())
        field(cell) = x_[cell.globalIndex()];
}

void EigenSparseMatrixSolver::mapSolution(VectorFiniteVolumeField &field)
{
    const Size nActiveCells = field.grid.nActiveCells();
    for(const Cell& cell: field.grid.activeCells())
    {
        Vector2D& vec = field(cell);
        vec.x = x_[cell.globalIndex()];
        vec.y = x_[cell.globalIndex() + nActiveCells];
    }
}

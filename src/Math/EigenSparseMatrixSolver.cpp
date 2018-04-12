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
    triplets.reserve(5 * coeffs.size());

    for (int i = 0, end = coeffs.size(); i < end; ++i)
        for (const auto &entry: coeffs[i])
            triplets.push_back(Triplet(i, entry.first, entry.second));

    mat_.setFromTriplets(triplets.begin(), triplets.end());
    mat_.makeCompressed();
}

void EigenSparseMatrixSolver::setGuess(const Vector &x0)
{
    for (int i = 0, end = x0.size(); i < end; ++i)
        x_(i) = x0(i);
}

void EigenSparseMatrixSolver::setRhs(const Vector &rhs)
{
    for (int i = 0, end = rhs.size(); i < end; ++i)
        rhs_(i) = rhs(i);
}

Scalar EigenSparseMatrixSolver::solve()
{
    solver_.compute(mat_);
    x_ = solver_.solve(rhs_);
    return 0.;
}

Scalar EigenSparseMatrixSolver::solve(const Vector &x0)
{
    return solve();
}
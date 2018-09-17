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

void EigenSparseMatrixSolver::setRank(int rowRank, int colRank)
{
    mat_.resize(rowRank, colRank);
    x_.resize(rowRank);
    rhs_.resize(colRank);
}

void EigenSparseMatrixSolver::set(const CoefficientList &coeffs)
{
    triplets_.clear();
    triplets_.reserve(5 * coeffs.size());

    for (int i = 0, end = coeffs.size(); i < end; ++i)
        for (const auto &entry: coeffs[i])
            triplets_.push_back(Triplet(i, entry.first, entry.second));

    mat_.setFromTriplets(triplets_.begin(), triplets_.end());
    mat_.makeCompressed();
}

void EigenSparseMatrixSolver::set(const std::vector<Index> &rowPtr, const std::vector<Index> &colInds, const std::vector<Scalar> &vals)
{
    triplets_.clear();
    triplets_.reserve(rowPtr.back());

    for(auto row = 0; row < rowPtr.size() - 1; ++row)
        for(auto j = rowPtr[row]; j < rowPtr[row + 1]; ++j)
            if(colInds[j] >= 0)
                triplets_.push_back(Triplet(row, colInds[j], vals[j]));

    mat_.setFromTriplets(triplets_.begin(), triplets_.end());
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

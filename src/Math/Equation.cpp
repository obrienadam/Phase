#include "Equation.h"

Equation::Equation(Size rank, Size nnz)
    :
      coeffs_(rank),
      rhs_(rank, 0.)
{
    std::for_each(coeffs_.begin(), coeffs_.end(), [nnz](SparseMatrixSolver::Row &row)
    { row.reserve(nnz); });
}

Equation &Equation::operator=(const Equation &eqn)
{
    if(this != &eqn)
    {
        solver_ = eqn.solver_ ? eqn.solver_: solver_;
        coeffs_ = eqn.coeffs_;
        rhs_ = eqn.rhs_;
    }

    return *this;
}

Equation &Equation::operator=(Equation &&eqn)
{
    solver_ = eqn.solver_ ? eqn.solver_: solver_;
    coeffs_ = std::move(eqn.coeffs_);
    rhs_ = std::move(eqn.rhs_);

    return *this;
}

void Equation::setRank(Size rank)
{
    coeffs_.resize(rank);
    rhs_.resize(rank, 0.);

    if(solver_)
        solver_->setRank(rank);
}

void Equation::setRank(Size nRows, Size nCols)
{
    coeffs_.resize(nRows);
    rhs_.resize(nRows);

    if(solver_)
        solver_->setRank(nRows, nCols);
}

void Equation::addRow(Size nnz)
{
    coeffs_.push_back(SparseMatrixSolver::Row());
    coeffs_.back().reserve(nnz);
    rhs_.resize(rhs_.size() + 1, 0.);
}

Scalar Equation::coeff(Index localRow, Index globalCol) const
{
    for (const SparseMatrixSolver::Entry &entry: coeffs_[localRow])
        if (entry.first == globalCol)
            return entry.second;

    return 0.;
}

Scalar &Equation::coeffRef(Index localRow, Index globalCol)
{
    for (SparseMatrixSolver::Entry &entry: coeffs_[localRow])
        if (entry.first == globalCol)
            return entry.second;
}

void Equation::setCoeff(Index localRow, Index globalCol, Scalar val)
{
    for (SparseMatrixSolver::Entry &entry: coeffs_[localRow])
        if (entry.first == globalCol)
        {
            entry.second = val;
            return;
        }

    coeffs_[localRow].push_back(SparseMatrixSolver::Entry(globalCol, val));
}

void Equation::addCoeff(Index localRow, Index globalCol, Scalar val)
{
    for (SparseMatrixSolver::Entry &entry: coeffs_[localRow])
        if (entry.first == globalCol)
        {
            entry.second += val;
            return;
        }

    coeffs_[localRow].push_back(SparseMatrixSolver::Entry(globalCol, val));
}

Scalar Equation::solve()
{
    solver_->setRank(coeffs_.size());
    solver_->set(coeffs_);
    solver_->setRhs(-rhs_);
    solver_->solve();
    return solver_->error();
}

Scalar Equation::solveLeastSquares()
{
    solver_->set(coeffs_);
    solver_->setRhs(-rhs_);
    solver_->solveLeastSquares();
    return solver_->error();
}

void Equation::clear()
{
    for (auto &row: coeffs_)
        row.clear();

    rhs_.zero();
}

Equation &Equation::operator+=(const Equation &rhs)
{
    for (int i = 0; i < coeffs_.size(); ++i)
        for (const auto &entry: rhs.coeffs_[i])
            addCoeff(i, entry.first, entry.second);

    rhs_ += rhs.rhs_;

    return *this;
}

Equation &Equation::operator-=(const Equation &rhs)
{
    for (int i = 0; i < coeffs_.size(); ++i)
        for (const auto &entry: rhs.coeffs_[i])
            addCoeff(i, entry.first, -entry.second);

    rhs_ -= rhs.rhs_;

    return *this;
}

Equation &Equation::operator+=(const Vector &rhs)
{
    rhs_ += rhs;
    return *this;
}

Equation &Equation::operator-=(const Vector &rhs)
{
    rhs_ -= rhs;
    return *this;
}

Equation &Equation::operator*=(Scalar rhs)
{
    for (int i = 0; i < coeffs_.size(); ++i)
        for (auto &entry: coeffs_[i])
            entry.second *= rhs;

    rhs_ *= rhs;

    return *this;
}

Equation &Equation::operator==(Scalar rhs)
{
    if (rhs != 0.)
        rhs_ -= rhs;

    return *this;
}

Equation &Equation::operator==(const Equation &rhs)
{
    for (int i = 0; i < coeffs_.size(); ++i)
        for (const auto &entry: rhs.coeffs_[i])
            addCoeff(i, entry.first, -entry.second);

    rhs_ -= rhs.rhs_;

    return *this;
}

Equation &Equation::operator==(const Vector &rhs)
{
    rhs_ -= rhs;
    return *this;
}

//- Protected methods

//- External functions

Equation operator+(Equation lhs, const Equation &rhs)
{
    lhs += rhs;
    return lhs;
}

Equation operator+(Equation lhs, const Vector &rhs)
{
    lhs += rhs;
    return lhs;
}

Equation operator-(Equation lhs, const Equation &rhs)
{
    lhs -= rhs;
    return lhs;
}

Equation operator-(Equation lhs, const Vector &rhs)
{
    lhs -= rhs;
    return lhs;
}

Equation operator*(Equation lhs, Scalar rhs)
{
    lhs *= rhs;
    return lhs;
}

Equation operator*(Scalar lhs, Equation rhs)
{
    return rhs * lhs;
}

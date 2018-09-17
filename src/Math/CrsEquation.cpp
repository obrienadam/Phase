#include <numeric>

#include "CrsEquation.h"

CrsEquation::CrsEquation(Size nRows, Size nnz)
    :
      rowPtr_(nRows + 1, nnz),
      colInd_(nRows * nnz, -1),
      vals_(nRows * nnz),
      rhs_(nRows, 0.)
{
    rowPtr_[0] = 0;
    std::partial_sum(rowPtr_.begin(), rowPtr_.end(), rowPtr_.begin());
}

CrsEquation &CrsEquation::operator=(const CrsEquation &eqn)
{
    if(this != &eqn)
    {
        solver_ = eqn.solver_ ? eqn.solver_: solver_;
        rowPtr_ = eqn.rowPtr_;
        colInd_ = eqn.colInd_;
        vals_ = eqn.vals_;
        rhs_ = eqn.rhs_;
    }

    return *this;
}

void CrsEquation::setRank(Size rank)
{
    rowPtr_.resize(rank + 1, rowPtr().back());
    rhs_.resize(rank, 0.);

    if(solver_)
        solver_->setRank(rank);
}

void CrsEquation::setRank(Size nRows, Size nCols)
{
    rowPtr_.resize(nRows + 1, rowPtr().back());
    rhs_.resize(nRows, 0.);

    if(solver_)
        solver_->setRank(nRows, nCols);
}

void CrsEquation::clear()
{
    rowPtr_.resize(1);
    colInd_.clear();
    vals_.clear();
    rhs_.clear();
}

void CrsEquation::addRow(Size nnz)
{
    rowPtr_.push_back(rowPtr().back() + nnz);
    colInd_.resize(colInd_.size() + nnz, -1);
    vals_.resize(vals_.size() + nnz);
    rhs_.resize(rhs_.size() + 1, 0.);
}

void CrsEquation::addRows(Size nRows, Size nnz)
{
    rowPtr_.resize(rowPtr_.size() + nRows, nnz);
    std::partial_sum(rowPtr_.end() - nRows - 1, rowPtr_.end(), rowPtr_.end() - nRows - 1);

    colInd_.resize(colInd_.size() + nnz * nRows, -1);
    vals_.resize(vals_.size() + nnz * nRows);
    rhs_.resize(rhs_.size() + nRows, 0.);
}

CrsEquation &CrsEquation::operator=(CrsEquation &&eqn)
{
    solver_ = eqn.solver_ ? eqn.solver_: solver_;
    rowPtr_ = std::move(eqn.rowPtr_);
    colInd_ = std::move(eqn.colInd_);
    vals_ = std::move(eqn.vals_);
    rhs_ = std::move(eqn.rhs_);

    return *this;
}

void CrsEquation::addCoeff(Index localRow, Index globalCol, Scalar val)
{
    for(auto j = rowPtr_[localRow]; j < rowPtr_[localRow + 1]; ++j)
    {
        if(colInd_[j] == globalCol)
        {
            vals_[j] += val;
            return;
        }
        else if(colInd_[j] == -1)
        {
            colInd_[j] = globalCol;
            vals_[j] = val;
            return;
        }
    }

    //- Did not find a suitable place, must resize sparse structure (potentially slow due to copies)
    colInd_.insert(colInd_.begin() + rowPtr_[localRow + 1], globalCol);
    vals_.insert(vals_.begin() + rowPtr_[localRow + 1], val);
    std::transform(rowPtr_.begin() + localRow + 1,
                   rowPtr_.end(),
                   rowPtr_.begin() + localRow + 1, [](Index i) { return i + 1; });
}

void CrsEquation::setCoeff(Index localRow, Index globalCol, Scalar val)
{
    for(auto j = rowPtr_[localRow]; j < rowPtr_[localRow + 1]; ++j)
    {
        if(colInd_[j] == globalCol)
        {
            vals_[j] = val;
            return;
        }
        else if(colInd_[j] == -1)
        {
            colInd_[j] = globalCol;
            vals_[j] = val;
            return;
        }
    }

    //- Did not find a suitable place, must resize sparse structure (potentially slow due to copies)
    colInd_.insert(colInd_.begin() + rowPtr_[localRow + 1], globalCol);
    vals_.insert(vals_.begin() + rowPtr_[localRow + 1], val);
    std::transform(rowPtr_.begin() + localRow + 1,
                   rowPtr_.end(),
                   rowPtr_.begin() + localRow + 1, [](Index i) { return i + 1; });
}

Scalar CrsEquation::coeff(Index localRow, Index globalCol) const
{
    for(auto j = rowPtr_[localRow]; j < rowPtr_[localRow + 1]; ++j)
        if(colInd_[j] == globalCol)
            return vals_[j];

    return 0.;
}

Scalar CrsEquation::solve()
{
    solver_->setRank(rank());
    solver_->set(rowPtr_, colInd_, vals_);
    solver_->setRhs(-rhs_);
    solver_->solve();
    return solver_->error();
}

Scalar CrsEquation::solveLeastSquares()
{
    solver_->set(rowPtr_, colInd_, vals_);
    solver_->setRhs(-rhs_);
    solver_->solveLeastSquares();
    return solver_->error();
}

//- Operators
CrsEquation& CrsEquation::operator +=(const CrsEquation &rhs)
{
    for(auto row = 0; row < rowPtr_.size() - 1; ++row)
        for(auto j = rhs.rowPtr_[row]; j < rhs.rowPtr_[row + 1]; ++j)
            addCoeff(row, rhs.colInd_[j], rhs.vals_[j]);

    rhs_ += rhs.rhs_;
    return *this;
}

CrsEquation& CrsEquation::operator -=(const CrsEquation &rhs)
{
    for(auto row = 0; row < rowPtr_.size() - 1; ++row)
        for(auto j = rhs.rowPtr_[row]; j < rhs.rowPtr_[row + 1]; ++j)
            addCoeff(row, rhs.colInd_[j], -rhs.vals_[j]);

    rhs_ -= rhs.rhs_;
    return *this;
}

CrsEquation &CrsEquation::operator+=(const Vector &rhs)
{
    rhs_ += rhs;
    return *this;
}

CrsEquation &CrsEquation::operator-=(const Vector &rhs)
{
    rhs_ -= rhs;
    return *this;
}

CrsEquation &CrsEquation::operator*=(Scalar rhs)
{
    std::transform(vals_.begin(), vals_.end(), vals_.begin(), [rhs](Scalar val) { return rhs * val; });
    return *this;
}

CrsEquation &CrsEquation::operator==(Scalar rhs)
{
    if(rhs != 0.)
        rhs_ -= rhs;
    return *this;
}

CrsEquation &CrsEquation::operator==(const CrsEquation &rhs)
{
    return operator -=(rhs);
}

CrsEquation &CrsEquation::operator==(const Vector &rhs)
{
    return operator -=(rhs);
}

std::ostream &operator<<(std::ostream &os, const CrsEquation &eqn)
{
    for(auto row = 0; row < eqn.rowPtr().size() - 1; ++row)
        for(auto j = eqn.rowPtr()[row]; j < eqn.rowPtr()[row + 1]; ++j)
            os << "(" << row << "," << eqn.colInd()[j] << "," << eqn.vals()[j] << ")\n";
    return os;
}

CrsEquation operator +(CrsEquation lhs, const CrsEquation &rhs)
{
    lhs += rhs;
    return lhs;
}

CrsEquation operator -(CrsEquation lhs, const CrsEquation &rhs)
{
    lhs -= rhs;
    return lhs;
}

CrsEquation operator +(CrsEquation lhs, const Vector &rhs)
{
    lhs += rhs;
    return lhs;
}

CrsEquation operator -(CrsEquation lhs, const Vector &rhs)
{
    lhs -= rhs;
    return lhs;
}

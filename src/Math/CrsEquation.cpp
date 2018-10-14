#include <numeric>

#include "CrsEquation.h"

std::vector<Index> CrsEquation::tmpRowPtr_, CrsEquation::tmpColInd_;

std::vector<Scalar> CrsEquation::tmpVals_;

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

Size CrsEquation::expand(Size row, Size nnz)
{
    colInd_.insert(colInd_.begin() + rowPtr_[row + 1], nnz, -1);
    vals_.insert(vals_.begin() + rowPtr_[row + 1], nnz, 0.);

    std::transform(rowPtr_.begin() + row + 1,
                   rowPtr_.end(),
                   rowPtr_.begin() + row + 1,
                   [nnz](Index i) { return i + nnz; });
}

Size CrsEquation::expand(Size nnz)
{
    std::vector<Index> newCols;
    std::vector<Scalar> newVals;

    newCols.reserve(colInd_.size() + nnz * rank());
    newVals.reserve(vals_.size() + nnz * rank());

    for(auto row = 0; row < rank(); ++row)
    {
        newCols.insert(newCols.end(),
                       colInd_.begin() + rowPtr_[row],
                       colInd_.begin() + rowPtr_[row + 1]);

        newCols.insert(newCols.end(), nnz, -1);

        newVals.insert(newVals.end(),
                       vals_.begin() + rowPtr_[row],
                       vals_.begin() + rowPtr_[row + 1]);

        newVals.insert(newVals.end(), nnz, 0.);
    }

    colInd_ = std::move(newCols);
    vals_ = std::move(newVals);

    int i = 1;
    std::transform(rowPtr_.begin() + 1, rowPtr_.end(), rowPtr_.begin() + 1,
                   [nnz, &i](Index idx) { return idx + nnz * i++; });
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

void CrsEquation::scaleRow(Index localRow, Scalar val)
{
    std::for_each(vals_.begin() + rowPtr_[localRow],
                  vals_.begin() + rowPtr_[localRow + 1],
            [val](Scalar &coeff) { coeff *= val; }
    );

    rhs_(localRow) *= val;
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
    tmpRowPtr_.clear();
    tmpColInd_.clear();
    tmpVals_.clear();

    tmpRowPtr_.emplace_back(0);

    for(auto row = 0; row < rank(); ++row)
    {
        tmpRowPtr_.emplace_back(tmpRowPtr_.back());

        for(auto j = rowPtr_[row]; j < rowPtr_[row + 1]; ++j)
        {
            if(vals_[j] == 0. || colInd_[j] < 0)
                continue;

            tmpColInd_.emplace_back(colInd_[j]);
            tmpVals_.emplace_back(vals_[j]);
            ++tmpRowPtr_.back();
        }

        for(auto j = rhs.rowPtr_[row]; j < rhs.rowPtr_[row + 1]; ++j)
        {
            if(rhs.vals_[j] == 0. || rhs.colInd_[j] < 0)
                continue;

            auto first = tmpColInd_.begin() + tmpRowPtr_[row];
            auto last = tmpColInd_.begin() + tmpRowPtr_[row + 1];

            auto it = std::find(first, last, rhs.colInd_[j]);

            if(it != last)
                tmpVals_[it - tmpColInd_.begin()] += rhs.vals_[j];
            else
            {
                tmpColInd_.emplace_back(rhs.colInd_[j]);
                tmpVals_.emplace_back(rhs.vals_[j]);
                ++tmpRowPtr_.back();
            }
        }
    }

    rowPtr_ = tmpRowPtr_;
    colInd_ = tmpColInd_;
    vals_ = tmpVals_;

    rhs_ += rhs.rhs_;
    return *this;
}

CrsEquation& CrsEquation::operator -=(const CrsEquation &rhs)
{
    tmpRowPtr_.clear();
    tmpColInd_.clear();
    tmpVals_.clear();

    tmpRowPtr_.emplace_back(0);

    for(auto row = 0; row < rank(); ++row)
    {
        tmpRowPtr_.emplace_back(tmpRowPtr_.back());

        for(auto j = rowPtr_[row]; j < rowPtr_[row + 1]; ++j)
        {
            if(vals_[j] == 0. || colInd_[j] < 0)
                continue;

            tmpColInd_.emplace_back(colInd_[j]);
            tmpVals_.emplace_back(vals_[j]);
            ++tmpRowPtr_.back();
        }

        for(auto j = rhs.rowPtr_[row]; j < rhs.rowPtr_[row + 1]; ++j)
        {
            if(rhs.vals_[j] == 0. || rhs.colInd_[j] < 0)
                continue;

            auto first = tmpColInd_.begin() + tmpRowPtr_[row];
            auto last = tmpColInd_.begin() + tmpRowPtr_[row + 1];

            auto it = std::find(first, last, rhs.colInd_[j]);

            if(it != last)
                tmpVals_[it - tmpColInd_.begin()] -= rhs.vals_[j];
            else
            {
                tmpColInd_.emplace_back(rhs.colInd_[j]);
                tmpVals_.emplace_back(-rhs.vals_[j]);
                ++tmpRowPtr_.back();
            }
        }
    }

    rowPtr_ = tmpRowPtr_;
    colInd_ = tmpColInd_;
    vals_ = tmpVals_;

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
    std::for_each(vals_.begin(), vals_.end(), [rhs](Scalar &val) { val *= rhs; });
    rhs_ *= rhs;
    return *this;
}

CrsEquation &CrsEquation::operator/=(Scalar rhs)
{
    std::for_each(vals_.begin(), vals_.end(), [rhs](Scalar &val) { val /= rhs; });
    rhs_ /= rhs;
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

CrsEquation operator *(CrsEquation lhs, Scalar rhs)
{
    lhs *= rhs;
    return lhs;
}

CrsEquation operator /(CrsEquation lhs, Scalar rhs)
{
    lhs /= rhs;
    return lhs;
}

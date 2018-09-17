#include <numeric>

#include "CooEquation.h"

std::vector<SparseEntry> CooEquation::tmp_;

CooEquation::CooEquation(Size nRows, Size nnz)
    :
      rhs_(nRows, 0.),
      rank_(0)
{
    entries_.reserve(nRows * nnz);
}

CooEquation &CooEquation::operator=(const CooEquation &eqn)
{
    if(this != &eqn)
    {
        solver_ = eqn.solver_ ? eqn.solver_: solver_;
        entries_ = eqn.entries_;
        rhs_ = eqn.rhs_;
        rank_ = eqn.rank_;
    }

    return *this;
}

CooEquation &CooEquation::operator=(CooEquation &&eqn)
{
    solver_ = eqn.solver_ ? eqn.solver_: solver_;
    entries_ = std::move(eqn.entries_);
    rhs_ = std::move(eqn.rhs_);
    rank_ = eqn.rank_;

    return *this;
}

void CooEquation::setRank(Size rank)
{
    rank_ = rank;
    rhs_.resize(rank_, 0.);

    if(solver_)
        solver_->setRank(rank);
}

void CooEquation::setRank(Size nRows, Size nCols)
{
    rank_ = nRows;

    rhs_.resize(rank_, 0.);

    if(solver_)
        solver_->setRank(nRows, nCols);
}

void CooEquation::clear()
{
    entries_.clear();
    rhs_.clear();
    rank_ = 0;
}

void CooEquation::addRow(Size nnz)
{
    addRows(1, nnz);
}

void CooEquation::addRows(Size nRows, Size nnz)
{
    entries_.reserve(entries_.size() + nRows * nnz);
    rhs_.resize(rhs_.size() + nRows, 0.);
    rank_ += nRows;
}

void CooEquation::addCoeff(Index localRow, Index globalCol, Scalar val)
{
    rank_ = std::max(rank_, (Size)localRow + 1);
    entries_.emplace_back(localRow, globalCol, val);
}

void CooEquation::setCoeff(Index localRow, Index globalCol, Scalar val)
{
    rank_ = std::max(rank_, (Size)localRow + 1);
    entries_.erase(std::remove_if(entries_.begin(), entries_.end(),
                                  [localRow, globalCol](const SparseEntry &e){ return e.row == localRow && e.col == globalCol; }),
                   entries_.end());

    addCoeff(localRow, globalCol, val);
}

Scalar CooEquation::coeff(Index localRow, Index globalCol) const
{
    Scalar val = 0.;

    for(const auto &e: entries_)
        if(e.row == localRow && e.col == globalCol)
            val += e.val;

    return val;
}

Scalar CooEquation::solve()
{
    solver_->setRank(rank());
    solver_->set(entries_);
    solver_->setRhs(-rhs_);
    solver_->solve();
    return solver_->error();
}

Scalar CooEquation::solveLeastSquares()
{
    solver_->set(entries_);
    solver_->setRhs(-rhs_);
    solver_->solveLeastSquares();
    return solver_->error();
}

//- Operators
CooEquation &CooEquation::operator +=(const std::vector<SparseEntry> &entries)
{
    entries_.insert(entries_.end(), entries.begin(), entries.end());
    return *this;
}

CooEquation& CooEquation::operator +=(const CooEquation &rhs)
{
    rhs_ += rhs.rhs_;
    return operator +=(rhs.entries_);
}

CooEquation& CooEquation::operator -=(const CooEquation &rhs)
{
    rhs_ -= rhs.rhs_;
    tmp_.assign(rhs.entries_.begin(), rhs.entries_.end());
    std::for_each(tmp_.begin(), tmp_.end(), [](SparseEntry &e) { e.val = -e.val; });
    return operator +=(tmp_);
}

CooEquation &CooEquation::operator+=(const Vector &rhs)
{
    rhs_ += rhs;
    return *this;
}

CooEquation &CooEquation::operator-=(const Vector &rhs)
{
    rhs_ -= rhs;
    return *this;
}

CooEquation &CooEquation::operator*=(Scalar rhs)
{
    std::for_each(entries_.begin(), entries_.end(), [rhs](SparseEntry &e) { e.val *= rhs; });
    rhs_ *= rhs;
    return *this;
}

CooEquation &CooEquation::operator==(Scalar rhs)
{
    if(rhs != 0.)
        rhs_ -= rhs;
    return *this;
}

CooEquation &CooEquation::operator==(const CooEquation &rhs)
{
    return operator -=(rhs);
}

CooEquation &CooEquation::operator==(const Vector &rhs)
{
    return operator -=(rhs);
}

std::ostream &operator<<(std::ostream &os, const CooEquation &eqn)
{
    for(const auto &e: eqn.entries())
        os << "(" << e.row << "," << e.col << "," << e.val << ")\n";

    return os;
}

CooEquation operator +(CooEquation lhs, const CooEquation &rhs)
{
    lhs += rhs;
    return lhs;
}

CooEquation operator -(CooEquation lhs, const CooEquation &rhs)
{
    lhs -= rhs;
    return lhs;
}

CooEquation operator +(CooEquation lhs, const Vector &rhs)
{
    lhs += rhs;
    return lhs;
}

CooEquation operator -(CooEquation lhs, const Vector &rhs)
{
    lhs -= rhs;
    return lhs;
}

#include "Equation.h"

Equation::Equation(Size rank, Size nnz)
    :
      _coeffs(rank),
      _rhs(rank, 0.)
{
    std::for_each(_coeffs.begin(), _coeffs.end(), [nnz](SparseMatrixSolver::Row &row)
    { row.reserve(nnz); });
}

Equation &Equation::operator=(const Equation &eqn)
{
    if(this != &eqn)
    {
        _spSolver = eqn._spSolver ? eqn._spSolver: _spSolver;
        _coeffs = eqn._coeffs;
        _rhs = eqn._rhs;
    }

    return *this;
}

Equation &Equation::operator=(Equation &&eqn)
{
    _spSolver = eqn._spSolver ? eqn._spSolver: _spSolver;
    _coeffs = std::move(eqn._coeffs);
    _rhs = std::move(eqn._rhs);

    return *this;
}

Scalar Equation::coeff(Index localRow, Index globalCol) const
{
    for (const SparseMatrixSolver::Entry &entry: _coeffs[localRow])
        if (entry.first == globalCol)
            return entry.second;

    return 0.;
}

Scalar &Equation::coeffRef(Index localRow, Index globalCol)
{
    for (SparseMatrixSolver::Entry &entry: _coeffs[localRow])
        if (entry.first == globalCol)
            return entry.second;
}

void Equation::clear()
{
    for (auto &row: _coeffs)
        row.clear();

    _rhs.zero();
}

Equation &Equation::operator+=(const Equation &rhs)
{
    for (int i = 0; i < _coeffs.size(); ++i)
        for (const auto &entry: rhs._coeffs[i])
            addCoeff(i, entry.first, entry.second);

    _rhs += rhs._rhs;

    return *this;
}

Equation &Equation::operator-=(const Equation &rhs)
{
    for (int i = 0; i < _coeffs.size(); ++i)
        for (const auto &entry: rhs._coeffs[i])
            addCoeff(i, entry.first, -entry.second);

    _rhs -= rhs._rhs;

    return *this;
}

Equation &Equation::operator+=(const Vector &rhs)
{
    _rhs += rhs;
    return *this;
}

Equation &Equation::operator-=(const Vector &rhs)
{
    _rhs -= rhs;
    return *this;
}

Equation &Equation::operator*=(Scalar rhs)
{
    for (int i = 0; i < _coeffs.size(); ++i)
        for (auto &entry: _coeffs[i])
            entry.second *= rhs;

    _rhs *= rhs;

    return *this;
}

Equation &Equation::operator==(Scalar rhs)
{
    if (rhs != 0.)
        _rhs -= rhs;

    return *this;
}

Equation &Equation::operator==(const Equation &rhs)
{
    for (int i = 0; i < _coeffs.size(); ++i)
        for (const auto &entry: rhs._coeffs[i])
            addCoeff(i, entry.first, -entry.second);

    _rhs -= rhs._rhs;

    return *this;
}

Equation &Equation::operator==(const Vector &rhs)
{
    _rhs -= rhs;
    return *this;
}

//- Protected methods

void Equation::setCoeff(Index localRow, Index globalCol, Scalar val)
{
    for (SparseMatrixSolver::Entry &entry: _coeffs[localRow])
        if (entry.first == globalCol)
        {
            entry.second = val;
            return;
        }

    _coeffs[localRow].push_back(SparseMatrixSolver::Entry(globalCol, val));
}

void Equation::addCoeff(Index localRow, Index globalCol, Scalar val)
{
    for (SparseMatrixSolver::Entry &entry: _coeffs[localRow])
        if (entry.first == globalCol)
        {
            entry.second += val;
            return;
        }

    _coeffs[localRow].push_back(SparseMatrixSolver::Entry(globalCol, val));
}

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

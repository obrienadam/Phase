#include "Term.h"

Term::Term(const FiniteVolumeGrid2D &grid)
{
    coefficients_.reserve(grid.nActiveCells());
    sources_.resize(grid.nActiveCells());
}

Term& Term::operator +=(const Term& rhs)
{
    for(int i = 0, end = coefficients_.size(); i < end; ++i)
    {
        int row = coefficients_[i].row();
        int col = coefficients_[i].col();
        Scalar value = coefficients_[i].value() + rhs.coefficients_[i].value();

        coefficients_[i] = Triplet(row, col, value);
    }

    for(int i = 0, end = sources_.size(); i < end; ++i)
        sources_[i] += rhs.sources_[i];

    return *this;
}

Term& Term::operator -=(const Term& rhs)
{
    for(int i = 0, end = coefficients_.size(); i < end; ++i)
    {
        int row = coefficients_[i].row();
        int col = coefficients_[i].col();
        Scalar value = coefficients_[i].value() - rhs.coefficients_[i].value();

        coefficients_[i] = Triplet(row, col, value);
    }

    for(int i = 0, end = sources_.size(); i < end; ++i)
        sources_[i] -= rhs.sources_[i];

    return *this;
}

Term& Term::operator *=(Scalar rhs)
{
    for(int i = 0, end = coefficients_.size(); i < end; ++i)
    {
        int row = coefficients_[i].row();
        int col = coefficients_[i].col();
        Scalar value = coefficients_[i].value()*rhs;

        coefficients_[i] = Triplet(row, col, value);
    }

    for(int i = 0, end = sources_.size(); i < end; ++i)
        sources_[i] *= rhs;

    return *this;
}

Term& Term::operator /=(Scalar rhs)
{
    for(int i = 0, end = coefficients_.size(); i < end; ++i)
    {
        int row = coefficients_[i].row();
        int col = coefficients_[i].col();
        Scalar value = coefficients_[i].value()/rhs;

        coefficients_[i] = Triplet(row, col, value);
    }

    for(int i = 0, end = sources_.size(); i < end; ++i)
        sources_[i] /= rhs;

    return *this;
}

//- External functions
Term operator==(Term term, Scalar rhs)
{
    for(int i = 0, end = term.sources_.size(); i < end; ++i)
    {
        term.sources_[i] += rhs;
    }

    return term;
}

Term operator+(Term lhs, const Term& rhs)
{
    lhs += rhs;
    return lhs;
}

Term operator-(Term lhs, const Term& rhs)
{
    lhs -= rhs;
    return lhs;
}

Term operator*(Term lhs, Scalar rhs)
{
    lhs *= rhs;
    return lhs;
}

Term operator*(Scalar lhs, Term rhs)
{
    rhs *= lhs;
    return rhs;
}

Term operator/(Term lhs, Scalar rhs)
{
    lhs /= rhs;
    return lhs;
}

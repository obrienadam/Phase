#include "Term.h"

Term::Term(const FiniteVolumeGrid2D &grid)
{
    coefficients_.reserve(grid.nCells());
    sources_.reserve(grid.nCells());
}

Term& Term::operator +=(const Term& rhs)
{
    for(int i = 0, end = coefficients_.size(); i < end; ++i)
    {
        coefficients_[i] += rhs.coefficients_[i];
        sources_[i] += rhs.sources_[i];
    }

    return *this;
}

Term& Term::operator -=(const Term& rhs)
{
    for(int i = 0, end = coefficients_.size(); i < end; ++i)
    {
        coefficients_[i] -= rhs.coefficients_[i];
        sources_[i] -= rhs.sources_[i];
    }

    return *this;
}

Term& Term::operator *=(Scalar rhs)
{
    for(int i = 0, end = coefficients_.size(); i < end; ++i)
    {
        coefficients_[i] *= rhs;
        sources_[i] *= rhs;
    }

    return *this;
}

Term& Term::operator /=(Scalar rhs)
{
    for(int i = 0, end = coefficients_.size(); i < end; ++i)
    {
        coefficients_[i] /= rhs;
        sources_[i] /= rhs;
    }

    return *this;
}

//- Triplet utility class
Term::Triplet& Term::Triplet::operator +=(const Term::Triplet& other)
{
    m_value += other.m_value;
    return *this;
}

Term::Triplet& Term::Triplet::operator -=(const Term::Triplet& other)
{
    m_value -= other.m_value;
    return *this;
}

Term::Triplet& Term::Triplet::operator *=(Scalar rhs)
{
    m_value *= rhs;
    return *this;
}

Term::Triplet& Term::Triplet::operator /=(Scalar rhs)
{
    m_value /= rhs;
    return *this;
}

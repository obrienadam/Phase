#include <assert.h>

#include "Term.h"

Term::Term(const ScalarFiniteVolumeField &var)
{
    coefficients_.reserve(var.grid.nActiveCells());
    sources_.resize(var.grid.nActiveCells());
}

Term::Term(const VectorFiniteVolumeField& var)
{
    coefficients_.reserve(2*var.grid.nActiveCells());
    sources_.resize(2*var.grid.nActiveCells());
}

Term& Term::operator +=(const Term& rhs)
{
    assert(sources_.size() == rhs.sources_.size());

    coefficients_.reserve(coefficients_.size() + rhs.coefficients_.size());

    for(const Triplet& entry: rhs.coefficients_)
        coefficients_.push_back(entry);

    for(int i = 0, end = sources_.size(); i < end; ++i)
        sources_[i] += rhs.sources_[i];

    return *this;
}

Term& Term::operator -=(const Term& rhs)
{
    assert(sources_.size() == rhs.sources_.size());

    coefficients_.reserve(coefficients_.size() + rhs.coefficients_.size());

    for(const Triplet& entry: rhs.coefficients_)
        coefficients_.push_back(Triplet(entry.row(), entry.col(), -entry.value()));

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

Term& Term::operator*=(const ScalarFiniteVolumeField& rhs)
{
    const size_t nScalars = rhs.size();

    assert(sources_.size()%nScalars == 0);

    for(int i = 0, end = coefficients_.size(); i < end; ++i)
    {
        int row = coefficients_[i].row();
        int col = coefficients_[i].col();
        Scalar value = coefficients_[i].value()*rhs[i%nScalars];

        coefficients_[i] = Triplet(row, col, value);
    }

    for(int i = 0, end = sources_.size(); i < end; ++i)
        sources_[i] *= rhs[i%nScalars];

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

Term operator==(Term term, const Term& rhs)
{
    term -= rhs;
    return term;
}

Term operator==(Term term, const ScalarFiniteVolumeField& field)
{
    for(int i = 0, end = term.sources_.size(); i < end; ++i)
    {
        term.sources_[i] += field[i];
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

Term operator*(const ScalarFiniteVolumeField& field, Term rhs)
{
    rhs *= field;
    return rhs;
}

Term operator*(Term lhs, const ScalarFiniteVolumeField& field)
{
    lhs *= field;
    return lhs;
}

Term operator/(Term lhs, Scalar rhs)
{
    lhs /= rhs;
    return lhs;
}

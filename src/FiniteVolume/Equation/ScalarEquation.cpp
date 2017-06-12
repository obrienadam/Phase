#include "Equation.h"

template<>
Equation<Scalar>::Equation(ScalarFiniteVolumeField &field, const std::string &name)
        :
        name(name),
        field_(field),
        nActiveCells_(field.grid().nLocalActiveCells()),
        coeffs_(nActiveCells_),
        sources_(nActiveCells_)
{
    for (auto &coeff: coeffs_)
        coeff.reserve(5);
}

template<>
template<>
void Equation<Scalar>::set(const Cell &cell, const Cell &nb, Scalar val)
{
    setValue(cell.index(0), nb.index(1), val);
}

template<>
template<>
void Equation<Scalar>::add(const Cell &cell, const Cell &nb, Scalar val)
{
    addValue(cell.index(0), nb.index(1), val);
}

template<>
Scalar Equation<Scalar>::get(const Cell &cell, const Cell &nb)
{
    for (const auto &entry: coeffs_[cell.index(0)])
        if (entry.first == nb.index(1))
            return entry.second;

    return 0.;
}

template<>
void Equation<Scalar>::remove(const Cell& cell)
{
    coeffs_[cell.index(0)].clear();
    sources_[cell.index(0)] = 0.;
}

template<>
void Equation<Scalar>::addSource(const Cell &cell, Scalar val)
{
    sources_[cell.index(0)] += val;
}

template<>
void Equation<Scalar>::setSource(const Cell &cell, Scalar val)
{
    sources_[cell.index(0)] = val;
}

template<>
void Equation<Scalar>::relax(Scalar relaxationFactor)
{
    for (const Cell &cell: field_.grid().localActiveCells())
    {
        const Index row = cell.index(0);
        Scalar &coeff = coeffRef(row, row);

        coeff /= relaxationFactor;
        sources_(row) -= (1. - relaxationFactor) * coeff * field_(cell);
    }
}

template<>
Equation<Scalar> &Equation<Scalar>::operator+=(const ScalarFiniteVolumeField &rhs)
{
    for (const Cell &cell: rhs.grid().localActiveCells())
        sources_(cell.index(0)) += rhs(cell);

    return *this;
}

template<>
Equation<Scalar> &Equation<Scalar>::operator-=(const ScalarFiniteVolumeField &rhs)
{
    for (const Cell &cell: rhs.grid().localActiveCells())
        sources_(cell.index(0)) -= rhs(cell);

    return *this;
}

//- Private
template<>
Size Equation<Scalar>::getRank() const
{
    return nActiveCells_;
}

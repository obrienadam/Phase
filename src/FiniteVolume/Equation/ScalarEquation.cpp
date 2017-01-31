#include "Equation.h"

template<>
Equation<Scalar>::Equation(ScalarFiniteVolumeField& field, const std::string& name)
    :
      name(name),
      field_(field),
      nActiveCells_(field.grid.nLocalActiveCells()),
      coeffs_(nActiveCells_),
      boundaries_(nActiveCells_),
      sources_(nActiveCells_)
{
    for(auto& coeff: coeffs_)
        coeff.reserve(5);
}

template<>
template<>
void Equation<Scalar>::set(const Cell &cell, const Cell &nb, Scalar val)
{
    setValue(cell.localIndex(), nb.globalIndex(0), val);
}

template<>
template<>
void Equation<Scalar>::add(const Cell &cell, const Cell &nb, Scalar val)
{
    addValue(cell.localIndex(), nb.globalIndex(0), val);
}

template<>
Scalar Equation<Scalar>::get(const Cell& cell, const Cell& nb)
{
    for(const auto& entry: coeffs_[cell.localIndex()])
        if(entry.first == nb.globalIndex(0))
            return entry.second;

    return 0.;
}


template<>
void Equation<Scalar>::addBoundary(const Cell& cell, Scalar val)
{
    boundaries_[cell.localIndex()] += val;
}

template<>
void Equation<Scalar>::setBoundary(const Cell& cell, Scalar val)
{
    boundaries_[cell.localIndex()] = val;
}

template<>
void Equation<Scalar>::addSource(const Cell& cell, Scalar val)
{
    sources_[cell.localIndex()] += val;
}

template<>
void Equation<Scalar>::setSource(const Cell& cell, Scalar val)
{
    sources_[cell.localIndex()] = val;
}

template<>
void Equation<Scalar>::relax(Scalar relaxationFactor)
{
    for(const Cell& cell: field_.grid.localActiveCells())
    {
        const Index row = cell.localIndex();
        Scalar &coeff = coeffRef(row, row);

        coeff /= relaxationFactor;
        boundaries_(row) += (1. - relaxationFactor)*coeff*field_(cell);
    }
}

template<>
Equation<Scalar>& Equation<Scalar>::operator +=(const ScalarFiniteVolumeField& rhs)
{
    for(const Cell& cell: rhs.grid.localActiveCells())
        sources_(cell.localIndex()) += rhs(cell);

    return *this;
}

template<>
Equation<Scalar>& Equation<Scalar>::operator -=(const ScalarFiniteVolumeField& rhs)
{
    for(const Cell& cell: rhs.grid.localActiveCells())
        sources_(cell.localIndex()) -= rhs(cell);

    return *this;
}

//- Private
template<>
Size Equation<Scalar>::getRank() const
{
    return nActiveCells_;
}

#include "Equation.h"

template<>
Equation<Scalar>::Equation(ScalarFiniteVolumeField &field, const std::string &name)
        :
        name(name),
        field_(field),
        coeffs_(field.grid()->localCells().size()),
        sources_(field.grid()->localCells().size(), 0.)
{
    std::for_each(coeffs_.begin(), coeffs_.end(), [](Row &row)
    {
        row.reserve(5);
    });
}

template<>
void Equation<Scalar>::set(const Cell &cell, const Cell &nb, Scalar val)
{
    setValue(field_.indexMap()->local(cell, 0), field_.indexMap()->global(nb, 0), val);
}

template<>
void Equation<Scalar>::add(const Cell &cell, const Cell &nb, Scalar val)
{
    addValue(field_.indexMap()->local(cell, 0), field_.indexMap()->global(nb, 0), val);
}

template<>
Scalar Equation<Scalar>::get(const Cell &cell, const Cell &nb)
{
    Index col = field_.indexMap()->global(nb, 0);

    for (const auto &entry: coeffs_[field_.indexMap()->local(cell, 0)])
        if (entry.first == col)
            return entry.second;

    return 0.;
}

template<>
void Equation<Scalar>::remove(const Cell &cell)
{
    coeffs_[field_.indexMap()->local(cell, 0)].clear();
    sources_(field_.indexMap()->local(cell, 0)) = 0.;
}

template<>
void Equation<Scalar>::addSource(const Cell &cell, Scalar val)
{
    sources_(field_.indexMap()->local(cell, 0)) += val;
}

template<>
void Equation<Scalar>::setSource(const Cell &cell, Scalar val)
{
    sources_(field_.indexMap()->local(cell, 0)) = val;
}

template<>
void Equation<Scalar>::relax(Scalar relaxationFactor)
{
    for (const Cell &cell: field_.grid()->localCells())
    {
        Index row = field_.indexMap()->local(cell, 0);
        Scalar &coeff = coeffRef(row, row);

        coeff /= relaxationFactor;
        sources_(row) -= (1. - relaxationFactor) * coeff * field_(cell);
    }
}

template<>
Equation<Scalar> &Equation<Scalar>::operator+=(const ScalarFiniteVolumeField &rhs)
{
    for (const Cell &cell: rhs.grid()->localCells())
        sources_(field_.indexMap()->local(cell, 0)) += rhs(cell);

    return *this;
}

template<>
Equation<Scalar> &Equation<Scalar>::operator-=(const ScalarFiniteVolumeField &rhs)
{
    for (const Cell &cell: rhs.grid()->localCells())
        sources_(field_.indexMap()->local(cell, 0)) -= rhs(cell);

    return *this;
}

//- Private

template<>
void Equation<Scalar>::mapFromSparseSolver()
{
    const IndexMap &idxMap = *field_.indexMap();

    for(const Cell& cell: field_.grid()->localCells())
        field_(cell) = sparseSolver()->x(idxMap.local(cell, 0));
}

template<>
Size Equation<Scalar>::getRank() const
{
    return field_.grid()->localCells().size();
}

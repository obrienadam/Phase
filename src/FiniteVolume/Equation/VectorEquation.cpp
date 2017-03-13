#include "Equation.h"

template<>
Equation<Vector2D>::Equation(VectorFiniteVolumeField &field, const std::string &name)
        :
        name(name),
        field_(field),
        coeffs_(2 * nActiveCells_),
        boundaries_(2 * nActiveCells_),
        sources_(2 * nActiveCells_)
{
    for (auto &coeff: coeffs_)
        coeff.reserve(5);
}

template<>
template<>
void Equation<Vector2D>::set(const Cell &cell, const Cell &nb, Scalar val)
{
    setValue(cell.index(0),
             nb.index(2),
             val);

    setValue(cell.index(0) + nActiveCells_,
             nb.index(3),
             val);
}

template<>
template<>
void Equation<Vector2D>::add(const Cell &cell, const Cell &nb, Scalar val)
{
    addValue(cell.index(0),
             nb.index(2),
             val);

    addValue(cell.index(0) + nActiveCells_,
             nb.index(3),
             val);
}

template<>
template<>
void Equation<Vector2D>::set(const Cell &cell, const Cell &nb, const Vector2D &val)
{
    setValue(cell.index(0),
             nb.index(2),
             val.x);

    setValue(cell.index(0) + nActiveCells_,
             nb.index(3),
             val.y);
}

template<>
template<>
void Equation<Vector2D>::add(const Cell &cell, const Cell &nb, const Vector2D &val)
{
    addValue(cell.index(0),
             nb.index(2),
             val.x);

    addValue(cell.index(0) + nActiveCells_,
             nb.index(3),
             val.y);
}

template<>
Vector2D Equation<Vector2D>::get(const Cell &cell, const Cell &nb)
{
    int i = 0;
    for (const auto &entry: coeffs_[cell.index(0)])
    {
        if (entry.first == nb.index(2))
            return Vector2D(entry.second, coeffs_[cell.index(0) + nActiveCells_][i].second);
        ++i;
    }

    return Vector2D(0., 0.);
}

template<>
void Equation<Vector2D>::addBoundary(const Cell &cell, Vector2D val)
{
    boundaries_[cell.index(0)] += val.x;
    boundaries_[cell.index(0) + nActiveCells_] += val.y;
}

template<>
void Equation<Vector2D>::relax(Scalar relaxationFactor)
{
    for (const Cell &cell: field_.grid.localActiveCells())
    {
        Scalar &coeffX = coeffRef(cell.index(0), cell.index(2));
        Scalar &coeffY = coeffRef(cell.index(0) + nActiveCells_, cell.index(3));

        coeffX /= relaxationFactor;
        coeffY /= relaxationFactor;

        boundaries_(cell.index(0)) += (1. - relaxationFactor) * coeffX * field_(cell).x;
        boundaries_(cell.index(0) + nActiveCells_) += (1. - relaxationFactor) * coeffY * field_(cell).y;
    }
}

template<>
Equation<Vector2D> &Equation<Vector2D>::operator+=(const VectorFiniteVolumeField &rhs)
{
    for (const Cell &cell: rhs.grid.localActiveCells())
    {
        Index rowX = cell.index(0);
        Index rowY = rowX + nActiveCells_;

        sources_(rowX) += rhs(cell).x;
        sources_(rowY) += rhs(cell).y;
    }

    return *this;
}

template<>
Equation<Vector2D> &Equation<Vector2D>::operator-=(const VectorFiniteVolumeField &rhs)
{
    for (const Cell &cell: rhs.grid.localActiveCells())
    {
        Index rowX = cell.index(0);
        Index rowY = rowX + nActiveCells_;

        sources_(rowX) -= rhs(cell).x;
        sources_(rowY) -= rhs(cell).y;
    }

    return *this;
}

//- Private
template<>
Size Equation<Vector2D>::getRank() const
{
    return 2 * nActiveCells_;
}

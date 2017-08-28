#include "Equation.h"

template<>
Equation<Vector2D>::Equation(VectorFiniteVolumeField &field, const std::string &name)
        :
        name(name),
        field_(field),
        nLocalActiveCells_(field.grid().nLocalActiveCells()),
        nGlobalActiveCells_(field.grid().nActiveCellsGlobal()),
        coeffs_(2 * nLocalActiveCells_),
        sources_(2 * nLocalActiveCells_)
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

    setValue(cell.index(0) + nLocalActiveCells_,
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

    addValue(cell.index(0) + nLocalActiveCells_,
             nb.index(3),
             val);
}

template<>
template<>
void Equation<Vector2D>::add(const Cell &cell, const Cell &nb, Vector2D val)
{
    addValue(cell.index(0),
             nb.index(2),
             val.x);

    addValue(cell.index(0) + nLocalActiveCells_,
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

    addValue(cell.index(0) + nLocalActiveCells_,
             nb.index(3),
             val.y);
}

template<>
template<>
void Equation<Vector2D>::addCoupling(const Cell &cell, const Cell &nb, const Vector2D &val)
{
    addValue(cell.index(0),
             nb.index(3),
             val.y);

    addValue(cell.index(0) + nLocalActiveCells_,
             nb.index(2),
             val.x);
}

template<>
Vector2D Equation<Vector2D>::get(const Cell &cell, const Cell &nb)
{
    Vector2D u(0., 0.);

    for (const auto &entry: coeffs_[cell.index(0)])
    {
        if (entry.first == nb.index(2))
        {
            u.x = entry.second;
            break;
        }
    }

    for (const auto &entry: coeffs_[cell.index(0) + nLocalActiveCells_])
    {
        if (entry.first == nb.index(3))
        {
            u.y = entry.second;
            break;
        }
    }

    return u;
}

template<>
void Equation<Vector2D>::remove(const Cell &cell)
{
    coeffs_[cell.index(0)].clear();
    coeffs_[cell.index(0) + nLocalActiveCells_].clear();
    sources_[cell.index(0)] = 0.;
    sources_[cell.index(0) + nLocalActiveCells_] = 0.;
}

template<>
void Equation<Vector2D>::addSource(const Cell& cell, Vector2D u)
{
    sources_[cell.index(0)] += u.x;
    sources_[cell.index(0) + nLocalActiveCells_] += u.y;
}

template<>
void Equation<Vector2D>::setSource(const Cell& cell, Vector2D u)
{
    sources_[cell.index(0)] = u.x;
    sources_[cell.index(0) + nLocalActiveCells_] = u.y;
}

template<>
void Equation<Vector2D>::relax(Scalar relaxationFactor)
{
    nLocalActiveCells_ = field_.grid().nLocalActiveCells();

    for (const Cell &cell: field_.grid().localActiveCells())
    {
        Scalar &coeffX = coeffRef(cell.index(0), cell.index(2));
        Scalar &coeffY = coeffRef(cell.index(0) + nLocalActiveCells_, cell.index(3));

        coeffX /= relaxationFactor;
        coeffY /= relaxationFactor;

        sources_(cell.index(0)) -= (1. - relaxationFactor) * coeffX * field_(cell).x;
        sources_(cell.index(0) + nLocalActiveCells_) -= (1. - relaxationFactor) * coeffY * field_(cell).y;
    }
}

template<>
Equation<Vector2D> &Equation<Vector2D>::operator+=(const VectorFiniteVolumeField &rhs)
{
    for (const Cell &cell: rhs.grid().localActiveCells())
    {
        Index rowX = cell.index(0);
        Index rowY = rowX + nLocalActiveCells_;

        sources_(rowX) += rhs(cell).x;
        sources_(rowY) += rhs(cell).y;
    }

    return *this;
}

template<>
Equation<Vector2D> &Equation<Vector2D>::operator-=(const VectorFiniteVolumeField &rhs)
{
    for (const Cell &cell: rhs.grid().localActiveCells())
    {
        Index rowX = cell.index(0);
        Index rowY = rowX + nLocalActiveCells_;

        sources_(rowX) -= rhs(cell).x;
        sources_(rowY) -= rhs(cell).y;
    }

    return *this;
}

//- Private
template<>
Size Equation<Vector2D>::getRank() const
{
    return 2 * field_.grid().localActiveCells().size();
}

#include "Geometry/Tensor2D.h"

#include "Equation.h"

template<>
Equation<Vector2D>::Equation(VectorFiniteVolumeField &field, const std::string &name)
        :
        name(name),
        field_(field),
        coeffs_(2 * field.grid()->localCells().size()),
        sources_(2 * field.grid()->localCells().size(), 0.)
{
    std::for_each(coeffs_.begin(), coeffs_.end(), [](Row &row)
    {
        row.reserve(5);
    });
}

template<>
void Equation<Vector2D>::set(const Cell &cell, const Cell &nb, Scalar val)
{
    setValue(field_.indexMap()->local(cell, 0),
             field_.indexMap()->global(nb, 0),
             val);

    setValue(field_.indexMap()->local(cell, 1),
             field_.indexMap()->global(nb, 1),
             val);
}

template<>
void Equation<Vector2D>::add(const Cell &cell, const Cell &nb, Scalar val)
{
    addValue(field_.indexMap()->local(cell, 0),
             field_.indexMap()->global(nb, 0),
             val);

    addValue(field_.indexMap()->local(cell, 1),
             field_.indexMap()->global(nb, 1),
             val);
}

template<>
template<>
void Equation<Vector2D>::add(const Cell &cell, const Cell &nb, const Vector2D &val)
{
    addValue(field_.indexMap()->local(cell, 0),
             field_.indexMap()->global(nb, 0),
             val.x);

    addValue(field_.indexMap()->local(cell, 1),
             field_.indexMap()->global(nb, 1),
             val.y);
}

template<>
template<>
void Equation<Vector2D>::add(const Cell &cell, const Cell &nb, const Tensor2D &val)
{
    addValue(field_.indexMap()->local(cell, 0),
             field_.indexMap()->global(nb, 0),
             val.xx);

    addValue(field_.indexMap()->local(cell, 0),
             field_.indexMap()->global(nb, 1),
             val.xy);

    addValue(field_.indexMap()->local(cell, 1),
             field_.indexMap()->global(nb, 0),
             val.yx);

    addValue(field_.indexMap()->local(cell, 1),
             field_.indexMap()->global(nb, 1),
             val.yy);
}

template<>
void Equation<Vector2D>::addCoupling(const Cell &cell, const Cell &nb, const Vector2D &val)
{
    addValue(field_.indexMap()->local(cell, 0),
             field_.indexMap()->global(nb, 1),
             val.y);

    addValue(field_.indexMap()->local(cell, 1),
             field_.indexMap()->global(nb, 0),
             val.x);
}

template<>
Vector2D Equation<Vector2D>::get(const Cell &cell, const Cell &nb)
{
    Vector2D u(0., 0.);

    Index col = field_.indexMap()->global(nb, 0);

    for (const auto &entry: coeffs_[field_.indexMap()->local(cell, 0)])
        if (entry.first == col)
        {
            u.x = entry.second;
            break;
        }

    col = field_.indexMap()->global(nb, 1);

    for (const auto &entry: coeffs_[field_.indexMap()->local(cell, 1)])
        if (entry.first == col)
        {
            u.y = entry.second;
            break;
        }

    return u;
}

template<>
void Equation<Vector2D>::remove(const Cell &cell)
{
    coeffs_[field_.indexMap()->local(cell, 0)].clear();
    coeffs_[field_.indexMap()->local(cell, 1)].clear();
    sources_(field_.indexMap()->local(cell, 0)) = 0.;
    sources_(field_.indexMap()->local(cell, 1)) = 0.;
}

template<>
void Equation<Vector2D>::addSource(const Cell &cell, Scalar val)
{
    sources_(field_.indexMap()->local(cell, 0)) += val;
    sources_(field_.indexMap()->local(cell, 1)) += val;
}

template<>
void Equation<Vector2D>::setSource(const Cell &cell, Scalar val)
{
    sources_(field_.indexMap()->local(cell, 0)) = val;
    sources_(field_.indexMap()->local(cell, 1)) = val;
}

template<>
template<>
void Equation<Vector2D>::addSource(const Cell &cell, const Vector2D &u)
{
    sources_(field_.indexMap()->local(cell, 0)) += u.x;
    sources_(field_.indexMap()->local(cell, 1)) += u.y;
}

template<>
template<>
void Equation<Vector2D>::setSource(const Cell &cell, const Vector2D &u)
{
    sources_(field_.indexMap()->local(cell, 0)) = u.x;
    sources_(field_.indexMap()->local(cell, 1)) = u.y;
}

template<>
template<>
void Equation<Vector2D>::addSource(const Cell &cell, const Tensor2D &tau)
{
    sources_(field_.indexMap()->local(cell, 0)) += tau.xx + tau.xy;
    sources_(field_.indexMap()->local(cell, 1)) += tau.yx + tau.yy;
}

template<>
void Equation<Vector2D>::relax(Scalar relaxationFactor)
{
    for (const Cell &cell: field_.grid()->localCells())
    {
        Scalar &coeffX = coeffRef(field_.indexMap()->local(cell, 0), field_.indexMap()->global(cell, 0));
        Scalar &coeffY = coeffRef(field_.indexMap()->local(cell, 1), field_.indexMap()->global(cell, 1));

        coeffX /= relaxationFactor;
        coeffY /= relaxationFactor;

        sources_(field_.indexMap()->local(cell, 0)) -= (1. - relaxationFactor) * coeffX * field_(cell).x;
        sources_(field_.indexMap()->local(cell, 1)) -= (1. - relaxationFactor) * coeffY * field_(cell).y;
    }
}

template<>
Equation<Vector2D> &Equation<Vector2D>::operator+=(const VectorFiniteVolumeField &rhs)
{
    for (const Cell &cell: rhs.grid()->localCells())
    {
        sources_(field_.indexMap()->local(cell, 0)) += rhs(cell).x;
        sources_(field_.indexMap()->local(cell, 1)) += rhs(cell).y;
    }

    return *this;
}

template<>
Equation<Vector2D> &Equation<Vector2D>::operator-=(const VectorFiniteVolumeField &rhs)
{
    for (const Cell &cell: rhs.grid()->localCells())
    {
        sources_(field_.indexMap()->local(cell, 0)) -= rhs(cell).x;
        sources_(field_.indexMap()->local(cell, 1)) -= rhs(cell).y;
    }

    return *this;
}

//- Private
template<>
void Equation<Vector2D>::mapFromSparseSolver()
{
    const IndexMap &idxMap = *field_.indexMap();

    for(const Cell& cell: field_.grid()->localCells())
    {
        field_(cell).x = sparseSolver()->x(idxMap.local(cell, 0));
        field_(cell).y = sparseSolver()->x(idxMap.local(cell, 1));
    }
}

template<>
Size Equation<Vector2D>::getRank() const
{
    return 2 * field_.grid()->localCells().size();
}

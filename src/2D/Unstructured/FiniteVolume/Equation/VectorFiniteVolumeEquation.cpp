#include "Geometry/Tensor2D.h"

#include "FiniteVolumeEquation.h"

template<>
FiniteVolumeEquation<Vector2D>::FiniteVolumeEquation(VectorFiniteVolumeField &field, const std::string &name)
        :
        Equation(2 * field.grid()->localCells().size(), 5),
        name(name),
        field_(field)
{

}

template<>
void FiniteVolumeEquation<Vector2D>::set(const Cell &cell, const Cell &nb, Scalar val)
{
    setCoeff(field_.indexMap()->local(cell, 0),
             field_.indexMap()->global(nb, 0),
             val);

    setCoeff(field_.indexMap()->local(cell, 1),
             field_.indexMap()->global(nb, 1),
             val);
}

template<>
void FiniteVolumeEquation<Vector2D>::add(const Cell &cell, const Cell &nb, Scalar val)
{
    addCoeff(field_.indexMap()->local(cell, 0),
             field_.indexMap()->global(nb, 0),
             val);

    addCoeff(field_.indexMap()->local(cell, 1),
             field_.indexMap()->global(nb, 1),
             val);
}

template<>
template<>
void FiniteVolumeEquation<Vector2D>::add(const Cell &cell, const Cell &nb, const Vector2D &val)
{
    addCoeff(field_.indexMap()->local(cell, 0),
             field_.indexMap()->global(nb, 0),
             val.x);

    addCoeff(field_.indexMap()->local(cell, 1),
             field_.indexMap()->global(nb, 1),
             val.y);
}

template<>
template<>
void FiniteVolumeEquation<Vector2D>::add(const Cell &cell, const Cell &nb, const Tensor2D &val)
{
    addCoeff(field_.indexMap()->local(cell, 0),
             field_.indexMap()->global(nb, 0),
             val.xx);

    addCoeff(field_.indexMap()->local(cell, 0),
             field_.indexMap()->global(nb, 1),
             val.xy);

    addCoeff(field_.indexMap()->local(cell, 1),
             field_.indexMap()->global(nb, 0),
             val.yx);

    addCoeff(field_.indexMap()->local(cell, 1),
             field_.indexMap()->global(nb, 1),
             val.yy);
}

template<>
void FiniteVolumeEquation<Vector2D>::addCoupling(const Cell &cell, const Cell &nb, const Vector2D &val)
{
    addCoeff(field_.indexMap()->local(cell, 0),
             field_.indexMap()->global(nb, 1),
             val.y);

    addCoeff(field_.indexMap()->local(cell, 1),
             field_.indexMap()->global(nb, 0),
             val.x);
}

template<>
Vector2D FiniteVolumeEquation<Vector2D>::get(const Cell &cell, const Cell &nb)
{
    Index rowX = field_.indexMap()->local(cell, 0);
    Index rowY = field_.indexMap()->local(cell, 1);
    Index colX = field_.indexMap()->global(nb, 0);
    Index colY = field_.indexMap()->global(nb, 1);

    return Vector2D(coeff(rowX, colX), coeff(rowY, colY));
}

template<>
void FiniteVolumeEquation<Vector2D>::remove(const Cell &cell)
{
    _coeffs[field_.indexMap()->local(cell, 0)].clear();
    _coeffs[field_.indexMap()->local(cell, 1)].clear();
    _rhs(field_.indexMap()->local(cell, 0)) = 0.;
    _rhs(field_.indexMap()->local(cell, 1)) = 0.;
}

template<>
void FiniteVolumeEquation<Vector2D>::addSource(const Cell &cell, Scalar val)
{
    _rhs(field_.indexMap()->local(cell, 0)) += val;
    _rhs(field_.indexMap()->local(cell, 1)) += val;
}

template<>
void FiniteVolumeEquation<Vector2D>::setSource(const Cell &cell, Scalar val)
{
    _rhs(field_.indexMap()->local(cell, 0)) = val;
    _rhs(field_.indexMap()->local(cell, 1)) = val;
}

template<>
template<>
void FiniteVolumeEquation<Vector2D>::addSource(const Cell &cell, const Vector2D &u)
{
    _rhs(field_.indexMap()->local(cell, 0)) += u.x;
    _rhs(field_.indexMap()->local(cell, 1)) += u.y;
}

template<>
template<>
void FiniteVolumeEquation<Vector2D>::setSource(const Cell &cell, const Vector2D &u)
{
    _rhs(field_.indexMap()->local(cell, 0)) = u.x;
    _rhs(field_.indexMap()->local(cell, 1)) = u.y;
}

template<>
template<>
void FiniteVolumeEquation<Vector2D>::addSource(const Cell &cell, const Tensor2D &tau)
{
    _rhs(field_.indexMap()->local(cell, 0)) += tau.xx + tau.xy;
    _rhs(field_.indexMap()->local(cell, 1)) += tau.yx + tau.yy;
}

template<>
void FiniteVolumeEquation<Vector2D>::relax(Scalar relaxationFactor)
{
    for (const Cell &cell: field_.grid()->localCells())
    {
        Scalar &coeffX = coeffRef(field_.indexMap()->local(cell, 0), field_.indexMap()->global(cell, 0));
        Scalar &coeffY = coeffRef(field_.indexMap()->local(cell, 1), field_.indexMap()->global(cell, 1));

        coeffX /= relaxationFactor;
        coeffY /= relaxationFactor;

        _rhs(field_.indexMap()->local(cell, 0)) -= (1. - relaxationFactor) * coeffX * field_(cell).x;
        _rhs(field_.indexMap()->local(cell, 1)) -= (1. - relaxationFactor) * coeffY * field_(cell).y;
    }
}

//- Private
template<>
void FiniteVolumeEquation<Vector2D>::mapFromSparseSolver()
{
    const IndexMap &idxMap = *field_.indexMap();

    for (const Cell &cell: field_.grid()->localCells())
    {
        field_(cell).x = _spSolver->x(idxMap.local(cell, 0));
        field_(cell).y = _spSolver->x(idxMap.local(cell, 1));
    }
}

template<>
Size FiniteVolumeEquation<Vector2D>::getRank() const
{
    return 2 * field_.grid()->localCells().size();
}

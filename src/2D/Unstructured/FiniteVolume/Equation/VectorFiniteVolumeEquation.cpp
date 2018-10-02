#include "Geometry/Tensor2D.h"

#include "FiniteVolumeEquation.h"

template<>
FiniteVolumeEquation<Vector2D>::FiniteVolumeEquation(VectorFiniteVolumeField &field, const std::string &name, int nnz)
    :
      CrsEquation(2 * field.grid()->localCells().size(), nnz),
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
void FiniteVolumeEquation<Vector2D>::scale(const Cell &cell, Scalar val)
{
    CrsEquation::scaleRow(field_.indexMap()->local(cell, 0), val);
    CrsEquation::scaleRow(field_.indexMap()->local(cell, 1), val);
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

    //- Added this to avoid coupling whenever possible (i.e. cartesian domains)
    if(val.xy != 0.)
        addCoeff(field_.indexMap()->local(cell, 0),
                 field_.indexMap()->global(nb, 1),
                 val.xy);

    if(val.yx != 0.)
        addCoeff(field_.indexMap()->local(cell, 1),
                 field_.indexMap()->global(nb, 0),
                 val.yx);

    addCoeff(field_.indexMap()->local(cell, 1),
             field_.indexMap()->global(nb, 1),
             val.yy);
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
    Index rowX = field_.indexMap()->local(cell, 0);
    Index rowY = field_.indexMap()->local(cell, 1);

    std::fill(colInd_.begin() + rowPtr_[rowX],
              colInd_.begin() + rowPtr_[rowX + 1],
            -1);

    std::fill(colInd_.begin() + rowPtr_[rowY],
              colInd_.begin() + rowPtr_[rowY + 1],
            -1);

    rhs_(rowX) = 0.;
    rhs_(rowY) = 0.;
}

template<>
void FiniteVolumeEquation<Vector2D>::addSource(const Cell &cell, Scalar val)
{
    rhs_(field_.indexMap()->local(cell, 0)) += val;
    rhs_(field_.indexMap()->local(cell, 1)) += val;
}

template<>
void FiniteVolumeEquation<Vector2D>::setSource(const Cell &cell, Scalar val)
{
    rhs_(field_.indexMap()->local(cell, 0)) = val;
    rhs_(field_.indexMap()->local(cell, 1)) = val;
}

template<>
template<>
void FiniteVolumeEquation<Vector2D>::addSource(const Cell &cell, const Vector2D &u)
{
    rhs_(field_.indexMap()->local(cell, 0)) += u.x;
    rhs_(field_.indexMap()->local(cell, 1)) += u.y;
}

template<>
template<>
void FiniteVolumeEquation<Vector2D>::setSource(const Cell &cell, const Vector2D &u)
{
    rhs_(field_.indexMap()->local(cell, 0)) = u.x;
    rhs_(field_.indexMap()->local(cell, 1)) = u.y;
}

template<>
template<>
void FiniteVolumeEquation<Vector2D>::addSource(const Cell &cell, const Tensor2D &tau)
{
    rhs_(field_.indexMap()->local(cell, 0)) += tau.xx + tau.xy;
    rhs_(field_.indexMap()->local(cell, 1)) += tau.yx + tau.yy;
}

template<>
void FiniteVolumeEquation<Vector2D>::relax(Scalar relaxationFactor)
{
    throw Exception("FiniteVolumeEquation<Vector2D>", "relax", "method removed.");

    //    for (const Cell &cell: field_.grid()->localCells())
    //    {
    //        Scalar &coeffX = coeffRef(field_.indexMap()->local(cell, 0), field_.indexMap()->global(cell, 0));
    //        Scalar &coeffY = coeffRef(field_.indexMap()->local(cell, 1), field_.indexMap()->global(cell, 1));

    //        coeffX /= relaxationFactor;
    //        coeffY /= relaxationFactor;

    //        _rhs(field_.indexMap()->local(cell, 0)) -= (1. - relaxationFactor) * coeffX * field_(cell).x;
    //        _rhs(field_.indexMap()->local(cell, 1)) -= (1. - relaxationFactor) * coeffY * field_(cell).y;
    //    }
}

//- Private
template<>
void FiniteVolumeEquation<Vector2D>::mapFromSparseSolver()
{
    const IndexMap &idxMap = *field_.indexMap();

    for (const Cell &cell: field_.grid()->localCells())
    {
        field_(cell).x = solver_->x(idxMap.local(cell, 0));
        field_(cell).y = solver_->x(idxMap.local(cell, 1));
    }
}

template<>
Size FiniteVolumeEquation<Vector2D>::getRank() const
{
    return 2 * field_.grid()->localCells().size();
}

//- Functions
template<>
FiniteVolumeEquation<Vector2D> operator * (const ScalarFiniteVolumeField &lhs, FiniteVolumeEquation<Vector2D> rhs)
{
    for(const Cell &cell: lhs.cells())
        rhs.scale(cell, lhs(cell));

    return rhs;
}

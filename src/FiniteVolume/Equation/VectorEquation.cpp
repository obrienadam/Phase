#include "Equation.h"

template<>
Equation<Vector2D>::Equation(VectorFiniteVolumeField& field, const std::string& name)
    :
      name(name),
      field_(field),
      nActiveCells_(field.grid.nActiveCells()),
      coeffs_(2*nActiveCells_),
      boundaries_(2*nActiveCells_),
      sources_(2*nActiveCells_)
{
    for(auto& coeff: coeffs_)
        coeff.reserve(10);
}

template<> template<>
void Equation<Vector2D>::set(const Cell &cell, const Cell &nb, Scalar val)
{
    setValue(cell.localIndex(),
             nb.globalIndex(1),
             val);

    setValue(cell.localIndex() + nActiveCells_,
             nb.globalIndex(2),
             val);
}

template<> template<>
void Equation<Vector2D>::add(const Cell &cell, const Cell &nb, Scalar val)
{
    addValue(cell.localIndex(),
             nb.globalIndex(1),
             val);

    addValue(cell.localIndex() + nActiveCells_,
             nb.globalIndex(2),
             val);
}

template<> template<>
void Equation<Vector2D>::set(const Cell &cell, const Cell &nb, const Vector2D& val)
{
    setValue(cell.localIndex(),
             nb.globalIndex(1),
             val.x);

    setValue(cell.localIndex() + nActiveCells_,
             nb.globalIndex(2),
             val.y);
}

template<> template<>
void Equation<Vector2D>::add(const Cell &cell, const Cell &nb, const Vector2D& val)
{
    addValue(cell.localIndex(),
             nb.globalIndex(1),
             val.x);

    addValue(cell.localIndex() + nActiveCells_,
             nb.globalIndex(2),
             val.y);
}

template<>
Vector2D Equation<Vector2D>::get(const Cell& cell, const Cell& nb)
{
    int i = 0;
    for(const auto& entry: coeffs_[cell.localIndex()])
    {
        if(entry.first == nb.globalIndex(1))
            return Vector2D(entry.second, coeffs_[cell.localIndex() + nActiveCells_][i].second);
        ++i;
    }

    return Vector2D(0., 0.);
}

template<>
void Equation<Vector2D>::addBoundary(const Cell& cell, Vector2D val)
{
    boundaries_[cell.localIndex()] += val.x;
    boundaries_[cell.localIndex() + nActiveCells_] += val.y;
}

//void setBoundary(const Cell& cell, T val);
//void addSource(const Cell& cell, T val);
//void setSource(const Cell& cell, T val)

template<>
void Equation<Vector2D>::relax(Scalar relaxationFactor)
{
    for(const Cell& cell: field_.grid.activeCells())
    {
        const Index rowX = cell.localIndex();
        const Index rowY = rowX + nActiveCells_;

        Scalar &coeffX = coeffRef(rowX, rowX);
        Scalar &coeffY = coeffRef(rowY, rowY);

        coeffX /= relaxationFactor;
        coeffY /= relaxationFactor;

        boundaries_(rowX) += (1. - relaxationFactor)*coeffX*field_(cell).x;
        boundaries_(rowY) += (1. - relaxationFactor)*coeffY*field_(cell).y;
    }
}

template<>
Equation<Vector2D>& Equation<Vector2D>::operator +=(const VectorFiniteVolumeField& rhs)
{
    for(const Cell& cell: rhs.grid.activeCells())
    {
        Index rowX = cell.localIndex();
        Index rowY = rowX + nActiveCells_;

        sources_(rowX) += rhs(cell).x;
        sources_(rowY) += rhs(cell).y;
    }

    return *this;
}

template<>
Equation<Vector2D>& Equation<Vector2D>::operator -=(const VectorFiniteVolumeField& rhs)
{
    for(const Cell& cell: rhs.grid.activeCells())
    {
        Index rowX = cell.localIndex();
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
    return 2*nActiveCells_;
}

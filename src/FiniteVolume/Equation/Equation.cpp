#include "Equation.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

template<>
Equation<ScalarFiniteVolumeField>::Equation(ScalarFiniteVolumeField& field, const std::string& name)
    :
      name(name),
      field_(field),
      coeffs_(field.grid.nActiveCells()),
      boundaries_(field.grid.nActiveCells()),
      sources_(field.grid.nActiveCells())
{
    for(auto& coeff: coeffs_)
        coeff.reserve(5);
}

template<>
Equation<VectorFiniteVolumeField>::Equation(VectorFiniteVolumeField& field, const std::string& name)
    :
      name(name),
      field_(field),
      coeffs_(2*field.grid.nActiveCells()),
      boundaries_(2*field.grid.nActiveCells()),
      sources_(2*field.grid.nActiveCells())
{
    for(auto& coeff: coeffs_)
        coeff.reserve(10);
}

template<>
Equation<ScalarFiniteVolumeField>& Equation<ScalarFiniteVolumeField>::operator +=(const ScalarFiniteVolumeField& rhs)
{
    for(const Cell& cell: rhs.grid.fluidCells())
        sources_(cell.globalIndex()) += rhs(cell);

    return *this;
}

template<>
Equation<ScalarFiniteVolumeField>& Equation<ScalarFiniteVolumeField>::operator -=(const ScalarFiniteVolumeField& rhs)
{
    for(const Cell& cell: rhs.grid.fluidCells())
        sources_(cell.globalIndex()) -= rhs(cell);

    return *this;
}

template<>
Equation<VectorFiniteVolumeField>& Equation<VectorFiniteVolumeField>::operator +=(const VectorFiniteVolumeField& rhs)
{
    const Size nActiveCells = rhs.grid.nActiveCells();

    for(const Cell& cell: rhs.grid.fluidCells())
    {
        Index rowX = cell.globalIndex();
        Index rowY = rowX + nActiveCells;

        sources_(rowX) += rhs(cell).x;
        sources_(rowY) += rhs(cell).y;
    }

    return *this;
}

template<>
Equation<VectorFiniteVolumeField>& Equation<VectorFiniteVolumeField>::operator -=(const VectorFiniteVolumeField& rhs)
{
    const Size nActiveCells = rhs.grid.nActiveCells();

    for(const Cell& cell: rhs.grid.fluidCells())
    {
        Index rowX = cell.globalIndex();
        Index rowY = rowX + nActiveCells;

        sources_(rowX) -= rhs(cell).x;
        sources_(rowY) -= rhs(cell).y;
    }

    return *this;
}

template<>
void Equation<ScalarFiniteVolumeField>::relax(Scalar relaxationFactor)
{
    for(const Cell& cell: field_.grid.fluidCells()) // Should the whole matrix be relaxed??
    {
        const Index row = cell.globalIndex();
        Scalar &coeff = getRef(row, row);

        boundaries_(row) += (1. - relaxationFactor)*coeff*field_(cell);
        coeff /= relaxationFactor;
    }
}

template<>
void Equation<VectorFiniteVolumeField>::relax(Scalar relaxationFactor)
{
    const Size nActiveCells = field_.grid.nActiveCells();

    for(const Cell& cell: field_.grid.fluidCells()) // Should the whole matrix be relaxed??
    {
        const Index rowX = cell.globalIndex();
        const Index rowY = rowX + nActiveCells;

        Scalar &coeffX = getRef(rowX, rowX);
        Scalar &coeffY = getRef(rowY, rowY);

        boundaries_(rowX) += (1. - relaxationFactor)*coeffX*field_(cell).x;
        boundaries_(rowY) += (1. - relaxationFactor)*coeffY*field_(cell).y;

        coeffX /= relaxationFactor;
        coeffY /= relaxationFactor;
    }
}

//- External functions

#include "Equation.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

template<>
Equation<ScalarFiniteVolumeField>::Equation(ScalarFiniteVolumeField& field, const std::string& name)
    :
      name(name),
      spMat_(field.grid.nActiveCells(), field.grid.nActiveCells()),
      field_(field),
      coeffs_(field.grid.nActiveCells())
{
    boundaries_ = SparseVector::Zero(field.grid.nActiveCells());
    sources_ = SparseVector::Zero(field.grid.nActiveCells());
    spMat_.reserve(5*field.grid.nActiveCells());

    for(auto& coeff: coeffs_)
        coeff.reserve(5);
}

template<>
Equation<VectorFiniteVolumeField>::Equation(VectorFiniteVolumeField& field, const std::string& name)
    :
      name(name),
      spMat_(2*field.grid.nActiveCells(), 2*field.grid.nActiveCells()),
      field_(field),
      coeffs_(2*field.grid.nActiveCells())
{
    boundaries_ = SparseVector::Zero(2*field.grid.nActiveCells());
    sources_ = SparseVector::Zero(2*field.grid.nActiveCells());
    spMat_.reserve(10*field.grid.nActiveCells());

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

        getRef(row, row) /= relaxationFactor;
        boundaries_(row) += (1. - relaxationFactor)*spMat_.coeff(row, row)*field_(cell);
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

        getRef(rowX, rowX) /= relaxationFactor;
        boundaries_(rowX) += (1. - relaxationFactor)*spMat_.coeff(rowX, rowX)*field_(cell).x;

        getRef(rowY, rowY) /= relaxationFactor;
        boundaries_(rowY) += (1. - relaxationFactor)*spMat_.coeff(rowY, rowY)*field_(cell).y;
    }
}

//- External functions

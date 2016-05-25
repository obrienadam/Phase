#include "Equation.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

template<>
Equation<ScalarFiniteVolumeField>::Equation(ScalarFiniteVolumeField& field, const std::string& name, SparseMatrix::Preconditioner precon)
    :
      name(name),
      spMat_(field.grid.nActiveCells(), field.grid.nActiveCells(), 5),
      field_(field),
      precon_(precon)
{
    boundaries_ = SparseVector::Zero(field.grid.nActiveCells());
    sources_ = SparseVector::Zero(field.grid.nActiveCells());
}

template<>
Equation<VectorFiniteVolumeField>::Equation(VectorFiniteVolumeField& field, const std::string& name, SparseMatrix::Preconditioner precon)
    :
      name(name),
      spMat_(2*field.grid.nActiveCells(), 2*field.grid.nActiveCells(), 5),
      field_(field),
      precon_(precon)
{
    boundaries_ = SparseVector::Zero(2*field.grid.nActiveCells());
    sources_ = SparseVector::Zero(2*field.grid.nActiveCells());
}

template<>
Equation<ScalarFiniteVolumeField>& Equation<ScalarFiniteVolumeField>::operator +=(const ScalarFiniteVolumeField& rhs)
{
    for(const Cell& cell: rhs.grid.fluidCells())
        sources_(cell.globalIndex()) += rhs[cell.id()];

    return *this;
}

template<>
Equation<ScalarFiniteVolumeField>& Equation<ScalarFiniteVolumeField>::operator -=(const ScalarFiniteVolumeField& rhs)
{
    for(const Cell& cell: rhs.grid.fluidCells())
        sources_(cell.globalIndex()) -= rhs[cell.id()];

    return *this;
}

template<>
Equation<VectorFiniteVolumeField>& Equation<VectorFiniteVolumeField>::operator +=(const VectorFiniteVolumeField& rhs)
{
    const size_t nActiveCells = rhs.grid.nActiveCells();

    for(const Cell& cell: rhs.grid.fluidCells())
    {
        size_t rowX = cell.globalIndex();
        size_t rowY = rowX + nActiveCells;

        sources_(rowX) += rhs[cell.id()].x;
        sources_(rowY) += rhs[cell.id()].y;
    }

    return *this;
}

template<>
Equation<VectorFiniteVolumeField>& Equation<VectorFiniteVolumeField>::operator -=(const VectorFiniteVolumeField& rhs)
{
    const size_t nActiveCells = rhs.grid.nActiveCells();

    for(const Cell& cell: rhs.grid.fluidCells())
    {
        size_t rowX = cell.globalIndex();
        size_t rowY = rowX + nActiveCells;

        sources_(rowX) -= rhs[cell.id()].x;
        sources_(rowY) -= rhs[cell.id()].y;
    }

    return *this;
}

template<>
void Equation<ScalarFiniteVolumeField>::relax(Scalar relaxationFactor)
{
    for(const Cell& cell: field_.grid.fluidCells()) // Should the whole matrix be relaxed??
    {
        size_t idx = cell.globalIndex();

        spMat_.coeffRef(idx, idx) /= relaxationFactor;
        boundaries_(idx) += (1. - relaxationFactor)*spMat_.coeff(idx, idx)*field_[cell.id()];
    }
}

template<>
void Equation<VectorFiniteVolumeField>::relax(Scalar relaxationFactor)
{
    const size_t nActiveCells = field_.grid.nActiveCells();

    for(const Cell& cell: field_.grid.fluidCells()) // Should the whole matrix be relaxed??
    {
        size_t idxX = cell.globalIndex();
        size_t idxY = idxX + nActiveCells;

        spMat_.coeffRef(idxX, idxX) /= relaxationFactor;
        boundaries_(idxX) += (1. - relaxationFactor)*spMat_.coeff(idxX, idxX)*field_[cell.id()].x;

        spMat_.coeffRef(idxY, idxY) /= relaxationFactor;
        boundaries_(idxY) += (1. - relaxationFactor)*spMat_.coeff(idxY, idxY)*field_[cell.id()].y;
    }
}

//- External functions

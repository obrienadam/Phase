#include "Equation.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

template<>
Equation<ScalarFiniteVolumeField>::Equation(const FiniteVolumeGrid2D &grid, ScalarFiniteVolumeField& field)
    :
      grid_(grid),
      spMat_(grid.nActiveCells(), grid.nActiveCells(), 5),
      b_(grid.nActiveCells()),
      field_(field)
{

}

template<>
Equation<VectorFiniteVolumeField>::Equation(const FiniteVolumeGrid2D &grid, VectorFiniteVolumeField& field)
    :
      grid_(grid),
      spMat_(2*grid.nActiveCells(), 2*grid.nActiveCells(), 5),
      b_(2*grid.nActiveCells()),
      field_(field)
{

}

template<>
void Equation<ScalarFiniteVolumeField>::relax(Scalar relaxationFactor)
{
    for(const Cell& cell: field_.grid.cells)
    {
        size_t idx = cell.globalIndex();

        spMat_.coeffRef(idx, idx) /= relaxationFactor;
        b_(idx) += (1. - relaxationFactor)*spMat_.coeff(idx, idx)*field_[cell.id()];
    }
}

template<>
void Equation<VectorFiniteVolumeField>::relax(Scalar relaxationFactor)
{
    const size_t nActiveCells = field_.grid.nActiveCells();

    for(const Cell& cell: field_.grid.cells)
    {
        size_t idxX = cell.globalIndex();
        size_t idxY = idxX + nActiveCells;

        spMat_.coeffRef(idxX, idxX) /= relaxationFactor;
        b_(idxX) += (1. - relaxationFactor)*spMat_.coeff(idxX, idxX)*field_[cell.id()].x;

        spMat_.coeffRef(idxY, idxY) /= relaxationFactor;
        b_(idxY) += (1. - relaxationFactor)*spMat_.coeff(idxY, idxY)*field_[cell.id()].y;
    }
}

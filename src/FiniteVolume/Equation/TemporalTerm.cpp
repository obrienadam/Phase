#include "TemporalTerm.h"

TemporalTerm::TemporalTerm(const ScalarFiniteVolumeField &var, Scalar timeStep)
    :
      Term(var)
{
    for(const Cell& cell: var.grid.cells)
    {
        size_t row = cell.globalIndex();
        Scalar coeff = cell.volume()/timeStep;

        coefficients_.push_back(Triplet(row, row, coeff));
        sources_[row] = coeff*var[cell.id()];
    }
}

TemporalTerm::TemporalTerm(const VectorFiniteVolumeField &var, Scalar timeStep)
    :
      Term(var)
{
    const size_t nActiveCells = var.grid.nActiveCells();

    for(const Cell& cell: var.grid.cells)
    {
        size_t rowX = cell.globalIndex();
        size_t rowY = rowX + nActiveCells;

        Scalar coeff = cell.volume()/timeStep;

        coefficients_.push_back(Triplet(rowX, rowX, coeff));
        coefficients_.push_back(Triplet(rowY, rowY, coeff));

        sources_[rowX] = coeff*var[cell.id()].x;
        sources_[rowY] = coeff*var[cell.id()].y;
    }
}

//- External functions

TemporalTerm ddt(const ScalarFiniteVolumeField &var, Scalar timeStep)
{
    return TemporalTerm(var, timeStep);
}

TemporalTerm ddt(const VectorFiniteVolumeField &var, Scalar timeStep)
{
    return TemporalTerm(var, timeStep);
}

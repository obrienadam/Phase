#include "TimeDerivative.h"

namespace fv
{
Equation<ScalarFiniteVolumeField> ddt(ScalarFiniteVolumeField& field, Scalar timeStep)
{
    const Field<Scalar> &prevField = field.prev(0);

    Equation<ScalarFiniteVolumeField> eqn(field);

    for(const Cell& cell: field.grid.cellZone("fluid"))
    {
        Index row = cell.localIndex();

        eqn.add(row, row, cell.volume()/timeStep);
        eqn.boundaries()(row) += cell.volume()*prevField[cell.id()]/timeStep;
    }

    return eqn;
}

Equation<VectorFiniteVolumeField> ddt(VectorFiniteVolumeField& field, Scalar timeStep)
{
    const size_t nActiveCells = field.grid.nActiveCells();
    const FiniteVolumeField<Vector2D> &prevField = field.prev(0);

    Equation<VectorFiniteVolumeField> eqn(field);

    for(const Cell& cell: field.grid.cellZone("fluid"))
    {
        Index rowX = cell.localIndex();
        Index rowY = rowX + nActiveCells;

        eqn.add(rowX, rowX, cell.volume()/timeStep);
        eqn.add(rowY, rowY, cell.volume()/timeStep);

        eqn.boundaries()(rowX) += cell.volume()*prevField[cell.id()].x/timeStep;
        eqn.boundaries()(rowY) += cell.volume()*prevField[cell.id()].y/timeStep;
    }

    return eqn;
}

Equation<VectorFiniteVolumeField> ddt(const ScalarFiniteVolumeField& rho, VectorFiniteVolumeField& field, Scalar timeStep)
{
    const size_t nActiveCells = field.grid.nActiveCells();
    const FiniteVolumeField<Vector2D> &prevField = field.prev(0);
    const ScalarFiniteVolumeField &rho0 = rho.prev(0);

    Equation<VectorFiniteVolumeField> eqn(field);

    for(const Cell& cell: field.grid.cellZone("fluid"))
    {
        const Index rowX = cell.localIndex();
        const Index rowY = rowX + nActiveCells;

        eqn.add(rowX, rowX, rho(cell)*cell.volume()/timeStep);
        eqn.add(rowY, rowY, rho(cell)*cell.volume()/timeStep);

        eqn.boundaries()(rowX) += rho0(cell)*cell.volume()*prevField(cell).x/timeStep;
        eqn.boundaries()(rowY) += rho0(cell)*cell.volume()*prevField(cell).y/timeStep;
    }

    return eqn;
}

}

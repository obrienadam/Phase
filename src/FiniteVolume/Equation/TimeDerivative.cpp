#include "TimeDerivative.h"

namespace fv
{
Equation<ScalarFiniteVolumeField> ddt(ScalarFiniteVolumeField& field, Scalar timeStep)
{
    const size_t nActiveCells = field.grid.nActiveCells();
    const Field<Scalar> &prevField = field.prev(0);

    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(field);

    entries.reserve(nActiveCells);

    for(const Cell& cell: field.grid.fluidCells())
    {
        size_t row = cell.globalIndex();

        entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, row, cell.volume()/timeStep));
        eqn.boundaries()(row) += cell.volume()*prevField[cell.id()]/timeStep;
    }

    eqn.assemble(entries);
    return eqn;
}

Equation<VectorFiniteVolumeField> ddt(VectorFiniteVolumeField& field, Scalar timeStep)
{
    const size_t nActiveCells = field.grid.nActiveCells();
    const FiniteVolumeField<Vector2D> &prevField = field.prev(0);

    std::vector<Equation<VectorFiniteVolumeField>::Triplet> entries;
    Equation<VectorFiniteVolumeField> eqn(field);

    entries.reserve(2*nActiveCells);

    for(const Cell& cell: field.grid.fluidCells())
    {
        size_t rowX = cell.globalIndex();
        size_t rowY = rowX + nActiveCells;

        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, rowX, cell.volume()/timeStep));
        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, rowY, cell.volume()/timeStep));

        eqn.boundaries()(rowX) += cell.volume()*prevField[cell.id()].x/timeStep;
        eqn.boundaries()(rowY) += cell.volume()*prevField[cell.id()].y/timeStep;
    }

    eqn.assemble(entries);
    return eqn;
}

Equation<VectorFiniteVolumeField> ddt(const ScalarFiniteVolumeField& rho, VectorFiniteVolumeField& field, Scalar timeStep)
{
    const size_t nActiveCells = field.grid.nActiveCells();
    const FiniteVolumeField<Vector2D> &prevField = field.prev(0);
    const ScalarFiniteVolumeField &rho0 = rho.prev(0);

    std::vector<Equation<VectorFiniteVolumeField>::Triplet> entries;
    Equation<VectorFiniteVolumeField> eqn(field);

    entries.reserve(2*nActiveCells);

    for(const Cell& cell: field.grid.fluidCells())
    {
        const Index rowX = cell.globalIndex();
        const Index rowY = rowX + nActiveCells;

        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, rowX, rho(cell)*cell.volume()/timeStep));
        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, rowY, rho(cell)*cell.volume()/timeStep));

        eqn.boundaries()(rowX) += rho0(cell)*cell.volume()*prevField(cell).x/timeStep;
        eqn.boundaries()(rowY) += rho0(cell)*cell.volume()*prevField(cell).y/timeStep;
    }

    eqn.assemble(entries);
    return eqn;
}

}

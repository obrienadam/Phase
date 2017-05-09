#include "GhostCellImmersedBoundaryObject.h"

GhostCellImmersedBoundaryObject::GhostCellImmersedBoundaryObject(const std::string &name, Label id, FiniteVolumeGrid2D &grid)
    :
      ImmersedBoundaryObject(name, id, grid)
{

}

void GhostCellImmersedBoundaryObject::update(Scalar timeStep)
{
    //- Doesn't move, stationary
}

void GhostCellImmersedBoundaryObject::updateCells()
{
    CellZone &fluidCells = grid_.cellZone("fluid");

    ibCells_->clear();
    solidCells_->clear();

    switch (shapePtr_->type())
    {
        case Shape2D::CIRCLE:
            for (const Cell &cell: fluidCells.cellCentersWithin(
                    *(Circle *) shapePtr_.get())) //- The circle method is much more efficient
                cells_->moveToGroup(cell);
            break;
        case Shape2D::BOX:
            for (const Cell &cell: fluidCells.cellCentersWithin(
                    *(Box *) shapePtr_.get())) //- The box method is much more efficient
                cells_->moveToGroup(cell);
            break;
        default:
            for (const Cell &cell: fluidCells.cellCentersWithin(*shapePtr_))
                cells_->moveToGroup(cell);
    }

    auto isIbCell = [this](const Cell &cell) {
        for (const InteriorLink &nb: cell.neighbours())
            if (!isInIb(nb.cell().centroid()))
                return true;

        //for (const DiagonalCellLink &dg: cell.diagonals())
        //    if (!isInIb(dg.cell().centroid()))
        //        return true;

        return false;
    };

    for (const Cell &cell: *cells_)
        if (isIbCell(cell))
            ibCells_->push_back(cell);
        else
            solidCells_->push_back(cell);

    grid_.setCellsActive(*ibCells_);
    grid_.setCellsInactive(*solidCells_);
    constructStencils();
}

Equation<Scalar> GhostCellImmersedBoundaryObject::bcs(ScalarFiniteVolumeField &field) const
{
    Equation<Scalar> eqn(field);
    BoundaryType bType = boundaryType(field.name());
    Scalar bRefValue = getBoundaryRefValue<Scalar>(field.name());

    for (const GhostCellStencil &st: stencils_)
    {
        Scalar centralCoeff;
        std::vector<Scalar> coeffs = st.ipCoeffs();

        //- Boundary assembly
        switch (bType)
        {
        case FIXED:
            centralCoeff = 0.5;
            for (Scalar &coeff: coeffs)
                coeff *= 0.5;

            eqn.addSource(st.cell(), -bRefValue);
            break;

        case ImmersedBoundaryObject::NORMAL_GRADIENT:
            centralCoeff = -1.;
            break;

        default:
            throw Exception("GhostCellImmersedBoundaryObject", "bcs", "invalid boundary type.");
        }

        eqn.add(st.cell(), st.cell(), centralCoeff);

        int i = 0;
        for (const Cell &ipCell: st.ipCells())
            eqn.add(st.cell(), ipCell, coeffs[i++]);
    }

    return eqn;
}

Equation<Vector2D> GhostCellImmersedBoundaryObject::bcs(VectorFiniteVolumeField &field) const
{
    Equation<Vector2D> eqn(field);
    BoundaryType bType = boundaryType(field.name());
    Vector2D bRefValue = getBoundaryRefValue<Vector2D>(field.name());

    for (const GhostCellStencil &st: stencils_)
    {
        Scalar centralCoeff;
        std::vector<Scalar> coeffs = st.ipCoeffs();

        //- Boundary assembly
        switch (bType)
        {
        case FIXED:
            centralCoeff = 0.5;
            for (Scalar &coeff: coeffs)
                coeff *= 0.5;

            eqn.addSource(st.cell(), -bRefValue);
            break;

        case ImmersedBoundaryObject::NORMAL_GRADIENT:
            centralCoeff = -1.;
            break;

        default:
            throw Exception("GhostCellImmersedBoundaryObject", "bcs", "invalid boundary type.");
        }

        eqn.add(st.cell(), st.cell(), centralCoeff);

        int i = 0;
        for (const Cell &ipCell: st.ipCells())
            eqn.add(st.cell(), ipCell, coeffs[i++]);
    }

    return eqn;
}

//- Protected methods

void GhostCellImmersedBoundaryObject::constructStencils()
{
    stencils_.clear();
    for (const Cell &cell: *ibCells_)
        stencils_.push_back(GhostCellStencil(cell, *shapePtr_, grid_));
}

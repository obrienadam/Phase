#include "ForcingCellImmersedBoundaryObject.h"

ForcingCellImmersedBoundaryObject::ForcingCellImmersedBoundaryObject(const std::string &name, Label id, FiniteVolumeGrid2D &grid)
    :
      ImmersedBoundaryObject(name, id, grid)
{

}

void ForcingCellImmersedBoundaryObject::update(Scalar timeStep)
{
    //- Does nothing, stationary
}

void ForcingCellImmersedBoundaryObject::updateCells()
{
    if(updateCellsCalled_)
        throw Exception("ForcingImmersedBoundaryObject", "updateCells", "cannot call updateCells more than once.");

    CellZone &fluidCells = grid_.cellZone("fluid");

    ibCells_->clear();
    solidCells_->clear();

    //- Find the solid cells
    switch (shapePtr_->type())
    {
    case Shape2D::CIRCLE:
        for (const Cell &cell: fluidCells.cellCentersWithin(
                 *(Circle *) shapePtr_.get())) //- The circle method is much more efficient
        {
            cells_->moveToGroup(cell);
            solidCells_->push_back(cell);
        }
        break;
    case Shape2D::BOX:
        for (const Cell &cell: fluidCells.cellCentersWithin(
                 *(Box *) shapePtr_.get())) //- The box method is much more efficient
        {
            cells_->moveToGroup(cell);
            solidCells_->push_back(cell);
        }
        break;
    default:
        for (const Cell &cell: fluidCells.cellCentersWithin(*shapePtr_))
        {
            cells_->moveToGroup(cell);
            solidCells_->push_back(cell);
        }
    }

    for(const Cell& cell: *solidCells_)
        for(const InteriorLink& nb: cell.neighbours())
            if(fluidCells.isInGroup(nb.cell()))
            {
                cells_->moveToGroup(nb.cell());
                ibCells_->push_back(nb.cell());
            }

    grid_.setCellsActive(*ibCells_);
    grid_.setCellsInactive(*solidCells_);
    constructStencils();
    updateCellsCalled_ = true;
}

Equation<Scalar> ForcingCellImmersedBoundaryObject::bcs(ScalarFiniteVolumeField &field) const
{
    Equation<Scalar> eqn(field);
    BoundaryType bType = boundaryType(field.name());
    Scalar bRefValue = getBoundaryRefValue<Scalar>(field.name());

    for (const ForcingCellStencil &st: stencils_)
    {
        //- Boundary assembly
        switch (bType)
        {
        case FIXED:
            eqn.add(st.cell(), st.cell(), 1.);
            eqn.add(st.cell(), st.iCells()[0], -st.iCoeffs()[0]);
            eqn.add(st.cell(), st.iCells()[1], -st.iCoeffs()[1]);
            eqn.addSource(st.cell(), -bRefValue*st.bCoeff());
            break;

        case NORMAL_GRADIENT:
            eqn.add(st.cell(), st.cell(), 1.);
            eqn.add(st.cell(), st.iCells()[0], -st.nCoeffs()[0]);
            eqn.add(st.cell(), st.iCells()[1], -st.nCoeffs()[1]);
            eqn.addSource(st.cell(), -bRefValue*st.bnCoeff()); //- Just to remind me how this works
            break;

        default:
            throw Exception("ForcedCellImmersedBoundaryObject", "bcs", "invalid boundary type.");
        }
    }

    return eqn;
}

Equation<Vector2D> ForcingCellImmersedBoundaryObject::bcs(VectorFiniteVolumeField &field) const
{
    Equation<Vector2D> eqn(field);
    BoundaryType bType = boundaryType(field.name());
    Vector2D bRefValue = getBoundaryRefValue<Vector2D>(field.name());

    for (const ForcingCellStencil &st: stencils_)
    {
        //- Boundary assembly
        switch (bType)
        {
        case FIXED:
            eqn.add(st.cell(), st.cell(), 1.);
            eqn.add(st.cell(), st.iCells()[0], -st.iCoeffs()[0]);
            eqn.add(st.cell(), st.iCells()[1], -st.iCoeffs()[1]);
            eqn.addSource(st.cell(), -bRefValue*st.bCoeff());
            break;

        case NORMAL_GRADIENT:
            eqn.add(st.cell(), st.cell(), 1.);
            eqn.add(st.cell(), st.iCells()[0], -st.nCoeffs()[0]);
            eqn.add(st.cell(), st.iCells()[1], -st.nCoeffs()[1]);
            eqn.addSource(st.cell(), -bRefValue*st.bnCoeff()); //- Just to remind me how this works
            break;

        default:
            throw Exception("ForcedCellImmersedBoundaryObject", "bcs", "invalid boundary type.");
        }
    }

    return eqn;
}

void ForcingCellImmersedBoundaryObject::constructStencils()
{
    stencils_.clear();
    for (const Cell &cell: *ibCells_)
        stencils_.push_back(ForcingCellStencil(cell, *shapePtr_, grid_.cellZone("fluid")));
}

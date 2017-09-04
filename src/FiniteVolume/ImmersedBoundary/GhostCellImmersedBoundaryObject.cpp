#include "GhostCellImmersedBoundaryObject.h"
#include "ContactLineGhostCellStencil.h"

GhostCellImmersedBoundaryObject::GhostCellImmersedBoundaryObject(const std::string &name, Label id,
                                                                 FiniteVolumeGrid2D &grid)
        :
        ImmersedBoundaryObject(name, id, grid)
{

}

void GhostCellImmersedBoundaryObject::update(Scalar timeStep)
{
    if (motion_)
    {
        motion_->update(*this, timeStep);
        updateCells();
    }
}

void GhostCellImmersedBoundaryObject::updateCells()
{
    if (motion_)
    {
        fluid_->add(freshCells_);
        freshCells_.clear();

        for (const Cell &cell: cells_)
            if (!isInIb(cell.centroid()))
                freshCells_.add(cell);
    }

    fluid_->add(cells_);

    ibCells_.clear();
    solidCells_.clear();
    deadCells_.clear();

    switch (shapePtr_->type())
    {
        case Shape2D::CIRCLE:
            for (const Cell &cell: fluid_->itemsWithin(
                    *std::static_pointer_cast<Circle>(shapePtr_))) //- The circle method is much more efficient
                cells_.add(cell);
            break;
        case Shape2D::BOX:
            for (const Cell &cell: fluid_->itemsWithin(
                    *std::static_pointer_cast<Box>(shapePtr_))) //- The box method is much more efficient
                cells_.add(cell);
            break;
        default:
            for (const Cell &cell: fluid_->itemsWithin(*shapePtr_))
                cells_.add(cell);
    }

    auto isIbCell = [this](const Cell &cell) {
        if (!isInIb(cell.centroid()))
            return false;

        for (const InteriorLink &nb: cell.neighbours())
            if (!isInIb(nb.cell().centroid()))
                return true;

        for (const DiagonalCellLink &dg: cell.diagonals())
            if (!isInIb(dg.cell().centroid()))
                return true;

        return false;
    };

    for (const Cell &cell: cells_)
        if (isIbCell(cell))
            ibCells_.add(cell);
        else
            solidCells_.add(cell);

    constructStencils();
}

Equation<Scalar> GhostCellImmersedBoundaryObject::bcs(ScalarFiniteVolumeField &field) const
{
    Equation<Scalar> eqn(field);
    BoundaryType bType = boundaryType(field.name());
    Scalar bRefValue = getBoundaryRefValue<Scalar>(field.name());

    for (const GhostCellStencil &st: stencils_)
    {
        //- Boundary assembly
        switch (bType)
        {
            case FIXED:
                eqn.add(st.cell(), st.cells(), st.dirichletCoeffs());
                eqn.addSource(st.cell(), -bRefValue);
                break;

            case ImmersedBoundaryObject::NORMAL_GRADIENT:
                eqn.add(st.cell(), st.cells(), st.neumannCoeffs());
                eqn.addSource(st.cell(), -bRefValue);
                break;

            default:
                throw Exception("GhostCellImmersedBoundaryObject", "bcs", "invalid boundary type.");
        }
    }

    for (const Cell &cell: solidCells_)
    {
        eqn.add(cell, cell, 1.);
        eqn.addSource(cell, -bRefValue);
    }

    return eqn;
}

Equation<Vector2D> GhostCellImmersedBoundaryObject::bcs(VectorFiniteVolumeField &field) const
{
    Equation<Vector2D> eqn(field);
    BoundaryType bType = boundaryType(field.name());
    Scalar bRefValue = getBoundaryRefValue<Scalar>(field.name());

    for (const GhostCellStencil &st: stencils_)
    {
        //- Boundary assembly
        switch (bType)
        {
            case FIXED:
                eqn.add(st.cell(), st.cells(), st.dirichletCoeffs());
                eqn.addSource(st.cell(), -bRefValue);
                break;

            case ImmersedBoundaryObject::NORMAL_GRADIENT:
                eqn.add(st.cell(), st.cells(), st.neumannCoeffs());
                eqn.addSource(st.cell(), -bRefValue);
                break;

            default:
                throw Exception("GhostCellImmersedBoundaryObject", "bcs", "invalid boundary type.");
        }
    }

    for (const Cell &cell: solidCells_)
    {
        eqn.add(cell, cell, 1.);
        eqn.addSource(cell, -bRefValue);
    }

    return eqn;
}

Equation<Vector2D> GhostCellImmersedBoundaryObject::solidVelocity(VectorFiniteVolumeField &u) const
{
    Equation<Vector2D> eqn(u);

    switch (boundaryType(u.name()))
    {
        case FIXED:
            for (const GhostCellStencil &st: stencils_)
            {
                eqn.add(st.cell(), st.cells(), st.dirichletCoeffs());
                eqn.addSource(st.cell(), -velocity(st.boundaryPoint()));
            }
            break;
        case PARTIAL_SLIP:
            break;
    }

    for (const Cell &cell: solidCells_)
    {
        eqn.add(cell, cell, 1.);
        eqn.addSource(cell, -velocity(cell.centroid()));
    }

    return eqn;
}

Equation<Scalar> GhostCellImmersedBoundaryObject::contactLineBcs(ScalarFiniteVolumeField& gamma, Scalar theta) const
{
    Equation<Scalar> eqn(gamma);

    for(const GhostCellStencil& st: stencils_)
    {
        auto clStencil = ContactLineGhostCellStencil(st.cell(), *this, gamma, theta);
        eqn.add(clStencil.cell(), clStencil.cells(), clStencil.neumannCoeffs());
    }

    for (const Cell &cell: solidCells_)
        eqn.add(cell, cell, 1.);

    return eqn;
}

//- Protected methods

void GhostCellImmersedBoundaryObject::constructStencils()
{
    stencils_.clear();
    for (const Cell &cell: ibCells_)
        stencils_.push_back(GhostCellStencil(cell, *this, grid_));
}

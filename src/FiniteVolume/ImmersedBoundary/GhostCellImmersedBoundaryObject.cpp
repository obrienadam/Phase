#include "GhostCellImmersedBoundaryObject.h"
#include "ForcingCellStencil.h"

GhostCellImmersedBoundaryObject::GhostCellImmersedBoundaryObject(const std::string &name, Label id, FiniteVolumeGrid2D &grid)
        :
        ImmersedBoundaryObject(name, id, grid)
{

}

void GhostCellImmersedBoundaryObject::update(Scalar timeStep)
{
    if(motion_)
    {
        motion_->update(*this, timeStep);
        updateCells();
    }
}

void GhostCellImmersedBoundaryObject::updateCells()
{
    fluid_->add(freshCells_); //- should remove from cells
    freshCells_.clear();

    for(const Cell& cell: cells_) //- New fresh cells
        if(!isInIb(cell.centroid()))
            freshCells_.add(cell);

    fluid_->add(cells_);
    ibCells_.clear();
    solidCells_.clear();

    switch (shapePtr_->type())
    {
        case Shape2D::CIRCLE:
            for (const Cell &cell: fluid_->itemsWithin(*std::static_pointer_cast<Circle>(shapePtr_))) //- The circle method is much more efficient
                cells_.add(cell);
            break;
        case Shape2D::BOX:
            for (const Cell &cell: fluid_->itemsWithin(*std::static_pointer_cast<Box>(shapePtr_))) //- The box method is much more efficient
                cells_.add(cell);
            break;
        default:
            for (const Cell &cell: fluid_->itemsWithin(*shapePtr_))
                cells_.add(cell);
    }

    auto isIbCell = [this](const Cell &cell) {
        if(!isInIb(cell.centroid()))
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
        std::vector<Scalar> coeffs = st.ipCoeffs();

        //- Boundary assembly
        switch (bType)
        {
            case FIXED:
                eqn.add(st.cell(), st.cell(), 1.);
                for(int i = 0; i < 4; ++i)
                    eqn.add(st.cell(), st.ipCells()[i], st.ipCoeffs()[i]);

                eqn.addSource(st.cell(), -2.*bRefValue);
                break;

            case ImmersedBoundaryObject::NORMAL_GRADIENT:
                eqn.add(st.cell(), st.cell(), -1.);
                for(int i = 0; i < 4; ++i)
                    eqn.add(st.cell(), st.ipCells()[i], st.ipCoeffs()[i]);
                break;

            default:
                throw Exception("GhostCellImmersedBoundaryObject", "bcs", "invalid boundary type.");
        }
    }

    for(const Cell& cell: solidCells_)
        eqn.add(cell, cell, 1.);

    return eqn;
}

Equation<Vector2D> GhostCellImmersedBoundaryObject::bcs(VectorFiniteVolumeField &field) const
{
    Equation<Vector2D> eqn(field);
    BoundaryType bType = boundaryType(field.name());
    Scalar bRefValue = getBoundaryRefValue<Scalar>(field.name());

    for (const GhostCellStencil &st: stencils_)
    {
        std::vector<Scalar> coeffs = st.ipCoeffs();

        //- Boundary assembly
        switch (bType)
        {
            case FIXED:
                eqn.add(st.cell(), st.cell(), 1.);
                for(int i = 0; i < 4; ++i)
                    eqn.add(st.cell(), st.ipCells()[i], st.ipCoeffs()[i]);

                eqn.addSource(st.cell(), -2.*bRefValue);
                break;

            case ImmersedBoundaryObject::NORMAL_GRADIENT:
                eqn.add(st.cell(), st.cell(), -1.);
                for(int i = 0; i < 4; ++i)
                    eqn.add(st.cell(), st.ipCells()[i], st.ipCoeffs()[i]);
                break;

            default:
                throw Exception("GhostCellImmersedBoundaryObject", "bcs", "invalid boundary type.");
        }
    }

    for(const Cell& cell: solidCells_)
        eqn.add(cell, cell, 1.);

    return eqn;
}

Equation<Vector2D> GhostCellImmersedBoundaryObject::solidVelocity(VectorFiniteVolumeField& u) const
{
    Equation<Vector2D> eqn(u);

    for(const GhostCellStencil& st: stencils_)
    {
        eqn.add(st.cell(), st.cell(), 1.);
        for(int i = 0; i < 4; ++i)
            eqn.add(st.cell(), st.ipCells()[i], st.ipCoeffs()[i]);

        eqn.addSource(st.cell(), -2*velocity(st.boundaryPoint()));
    }

    for(const Cell& cell: solidCells_)
    {
        eqn.add(cell, cell, 1.);
        eqn.addSource(cell, -velocity(cell.centroid()));
    }

    return eqn;
}

//- Protected methods

void GhostCellImmersedBoundaryObject::constructStencils()
{
    stencils_.clear();
    for (const Cell &cell: ibCells_)
        stencils_.push_back(GhostCellStencil(cell, *shapePtr_, grid_));
}

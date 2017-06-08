#include <fstream>

#include "ForcingCellImmersedBoundaryObject.h"

ForcingCellImmersedBoundaryObject::ForcingCellImmersedBoundaryObject(const std::string &name, Label id, FiniteVolumeGrid2D &grid)
    :
      ImmersedBoundaryObject(name, id, grid)
{
    pseudoFluidPoints_ = CellZone("PseudoFluidPoints", zoneRegistry_);
}

void ForcingCellImmersedBoundaryObject::update(Scalar timeStep)
{
    if(motion_)
    {
        motion_->update(*this, timeStep);
        updateCells();
    }
}

void ForcingCellImmersedBoundaryObject::updateCells()
{
    ibCells_.clear();
    solidCells_.clear();
    pseudoFluidPoints_.clear();

    fluid_->add(cells_);

    switch (shapePtr_->type())
    {
    case Shape2D::CIRCLE:
        for (const Cell &cell: fluid_->itemsWithin(
                 *(Circle *) shapePtr_.get())) //- The circle method is much more efficient
            cells_.add(cell);
        break;
    case Shape2D::BOX:
        for (const Cell &cell: fluid_->itemsWithin(
                 *(Box *) shapePtr_.get())) //- The box method is much more efficient
            cells_.add(cell);
        break;
    default:
        for (const Cell &cell: fluid_->itemsWithin(*shapePtr_))
            cells_.add(cell);
    }

    solidCells_.add(cells_);

    for(const Cell& cell: cells_.items())
        for(const InteriorLink& nb: cell.neighbours())
            if(!shapePtr_->isInside(nb.cell().centroid()))
            {
                cells_.add(nb.cell());
                ibCells_.add(nb.cell());
            }

    for(const Cell& cell: solidCells_.items()) //- Must make copy since container will be modified
        for(const InteriorLink& nb: cell.neighbours())
            if(ibCells_.isInGroup(nb.cell()))
            {
                pseudoFluidPoints_.add(cell);
                break;
            }

    grid_.setCellsActive(*fluid_);
    grid_.setCellsActive(ibCells_);
    grid_.setCellsActive(pseudoFluidPoints_);
    grid_.setCellsInactive(solidCells_);
    constructStencils();
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
            eqn.add(st.cell(), st.nbCells()[0], -st.dirichletCellCoeffs()[0]);
            eqn.add(st.cell(), st.nbCells()[1], -st.dirichletCellCoeffs()[1]);
            eqn.addSource(st.cell(), -st.dirichletBoundaryCoeff()*bRefValue);
            break;

        case NORMAL_GRADIENT:
            eqn.add(st.cell(), st.cell(), 1.);
            eqn.add(st.cell(), st.nbCells()[0], -st.neumannCellCoeffs()[0]);
            eqn.add(st.cell(), st.nbCells()[1], -st.neumannCellCoeffs()[1]);
            eqn.addSource(st.cell(), -st.neumannBoundaryCoeff()*bRefValue); //- Just to remind me how this works
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

    for (const ForcingCellStencil &st: stencils_)
    {
        Vector2D vel = velocity(st.xc());

        //- Boundary assembly
        switch (bType)
        {
        case FIXED:
            eqn.add(st.cell(), st.cell(), 1.);
            eqn.add(st.cell(), st.nbCells()[0], -st.dirichletCellCoeffs()[0]);
            eqn.add(st.cell(), st.nbCells()[1], -st.dirichletCellCoeffs()[1]);
            eqn.addSource(st.cell(), -st.dirichletBoundaryCoeff()*vel);
            break;

        case NORMAL_GRADIENT:
            eqn.add(st.cell(), st.cell(), 1.);
            eqn.add(st.cell(), st.nbCells()[0], -st.neumannCellCoeffs()[0]);
            eqn.add(st.cell(), st.nbCells()[1], -st.neumannCellCoeffs()[1]);
            eqn.addSource(st.cell(), -st.neumannBoundaryCoeff()*vel); //- Just to remind me how this works
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
    for (const Cell &cell: ibCells_)
        stencils_.push_back(ForcingCellStencil(cell, *shapePtr_, *fluid_));

    for(const Cell &cell: pseudoFluidPoints_)
        stencils_.push_back(ForcingCellStencil(cell, *shapePtr_, ibCells_));

    std::ofstream fout("stencils.dat");

    for(const ForcingCellStencil& st: stencils_)
    {
        fout << "Geometry x=0 y=0 T=LINE C=BLACK LT=0.2 CS=GRID\n"
             << "1\n"
             << "4\n"
             << st.xc().x << " " << st.xc().y << "\n"
             << st.nbCells()[0].get().centroid().x << " " << st.nbCells()[0].get().centroid().y << "\n"
             << st.nbCells()[1].get().centroid().x << " " << st.nbCells()[1].get().centroid().y << "\n"
             << st.xc().x << " " << st.xc().y << "\n";
    }

    fout.close();
}

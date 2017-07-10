#include "StepImmersedBoundaryObject.h"

void StepImmersedBoundaryObject::updateCells()
{
    fluid_->add(cells_);

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

    grid_.setCellsActive(*fluid_);
    grid_.setCellsActive(cells_);
}

Equation<Scalar> StepImmersedBoundaryObject::bcs(ScalarFiniteVolumeField& field) const
{
    Equation<Scalar> eqn(field);

    auto bType = boundaryType(field.name());
    auto bRefValue = getBoundaryRefValue<Scalar>(field.name());

    switch(bType)
    {
        case FIXED:
            for(const Cell& cell: cells_)
            {
                eqn.add(cell, cell, 1.);
                eqn.addSource(cell, -bRefValue);
            }
            break;
    }

    return eqn;
}

Equation<Vector2D> StepImmersedBoundaryObject::bcs(VectorFiniteVolumeField& field) const
{
    Equation<Vector2D> eqn(field);

    auto bType = boundaryType(field.name());

    switch(bType)
    {
        case FIXED:
            for(const Cell& cell: cells_)
            {
                eqn.add(cell, cell, 1.);
                eqn.addSource(cell, -velocity(cell.centroid()));
            }
            break;
    }

    return eqn;
}
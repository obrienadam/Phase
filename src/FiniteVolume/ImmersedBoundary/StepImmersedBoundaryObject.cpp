#include "StepImmersedBoundaryObject.h"

void StepImmersedBoundaryObject::update(Scalar timeStep)
{
    if(motion_)
    {
        motion_->update(*this, timeStep);
        updateCells();
    }
}

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

    solidCells_.clear();
    solidCells_.add(cells_);

    ibCells_.clear();
    for(const Cell& cell: solidCells_)
        for(const InteriorLink& nb: cell.neighbours())
            if(!isInIb(nb.cell().centroid()))
                ibCells_.add(nb.cell()); //- Note: these should still be fluid cells
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

void StepImmersedBoundaryObject::computeForce(Scalar rho,
                                              Scalar mu,
                                              const VectorFiniteVolumeField &u,
                                              const ScalarFiniteVolumeField &p)
{
    typedef std::pair<Point2D, Scalar> PressurePoint;

    std::vector<PressurePoint> pPoints;

    for(const Cell& cell: *fluid_)
    {
        for(const InteriorLink& nb: cell.neighbours())
        {
            if(isInIb(nb.cell().centroid()))
            {
                const Cell &stCell = fluid_->nearestItems(2*cell.centroid() - nb.cell().centroid(), 1)[0];
                LineSegment2D lnb = intersectionLine(LineSegment2D(cell.centroid(), nb.cell().centroid()));

                Vector2D rvA = lnb.ptB() - cell.centroid();
                Vector2D rvB = nb.cell().centroid() - lnb.ptB();

                Scalar g = rvB.mag()/(rvA.mag() + rvB.mag());

                Scalar pb = g*p(cell) + (1. - g)*p(nb.cell());

                pPoints.push_back(PressurePoint(lnb.ptB(), pb));
            }
        }
    }

    std::sort(pPoints.begin(), pPoints.end(), [this](const PressurePoint& ptA, const PressurePoint& ptB) {
        Vector2D rVecA = ptA.first - shapePtr_->centroid();
        Vector2D rVecB = ptB.first - shapePtr_->centroid();
        Scalar thetaA = atan2(rVecA.y, rVecA.x);
        Scalar thetaB = atan2(rVecB.y, rVecB.x);
        return (thetaA < 0. ? thetaA + 2*M_PI: thetaA) < (thetaB < 0. ? thetaB + 2*M_PI: thetaB);
    });

    force_ = Vector2D(0., 0.);
    if(shapePtr_->type() == Shape2D::CIRCLE)
    {
        auto circ = std::static_pointer_cast<Circle>(shapePtr_);
        Scalar r = circ->radius();

        for(int i = 0, end = pPoints.size(); i < end; ++i)
        {
            const PressurePoint &ptA = pPoints[i];
            const PressurePoint &ptB = pPoints[(i + 1) % end];

            Vector2D rA = ptA.first - circ->centroid();
            Vector2D rB = ptB.first - circ->centroid();
            Scalar a = atan2(rA.y, rA.x);
            Scalar b = atan2(rB.y, rB.x);
            a = a < 0 ? a + 2 * M_PI : a;
            b = b < 0 ? b + 2 * M_PI : b;

            if(fabs(b - a) < 1e-14)
                continue;

            Scalar m = (ptB.second - ptA.second) / (b - a);
            Scalar c = ptA.second - m * a;
            Scalar px = c * r * (sin(b) - sin(a)) + m * r * (-a * sin(a) - cos(a) + b * sin(b) + cos(b));
            Scalar py = r * cos(a) * (a * m + c) - r * (m * (sin(a) - sin(b)) + cos(b) * (b * m + c));

            force_ -= Vector2D(px, py);
        }
    }
}
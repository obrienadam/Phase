#include <math.h>

#include "StepImmersedBoundaryObject.h"
#include "Matrix.h"

void StepImmersedBoundaryObject::updateCells()
{
    fluid_->add(cells_);

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

    solidCells_.clear();
    solidCells_.add(cells_);

    ibCells_.clear();
    for (const Cell &cell: solidCells_)
        for (const InteriorLink &nb: cell.neighbours())
            if (!isInIb(nb.cell().centroid()))
                ibCells_.add(nb.cell()); //- Note: these should still be fluid cells
}

Equation<Scalar> StepImmersedBoundaryObject::bcs(ScalarFiniteVolumeField &field) const
{
    Equation<Scalar> eqn(field);

    auto bType = boundaryType(field.name());
    auto bRefValue = getBoundaryRefValue<Scalar>(field.name());

    switch (bType)
    {
        case FIXED:
            for (const Cell &cell: cells_)
            {
                eqn.add(cell, cell, 1.);
                eqn.addSource(cell, -bRefValue);
            }
            break;
    }

    return eqn;
}

Equation<Vector2D> StepImmersedBoundaryObject::bcs(VectorFiniteVolumeField &field) const
{
    Equation<Vector2D> eqn(field);

    auto bType = boundaryType(field.name());

    switch (bType)
    {
        case FIXED:
            for (const Cell &cell: cells_)
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
                                              const ScalarFiniteVolumeField &p,
                                              const Vector2D &g)
{
    typedef std::pair<Point2D, Scalar> ScalarPoint;
    typedef std::pair<Point2D, Vector2D> VectorPoint;

    std::vector<Point2D> bPts;
    std::vector<Scalar> bP;

    for (const Cell &cell: ibCells_)
        for (const InteriorLink &nb: cell.neighbours())
            if (isInIb(nb.cell()))
            {
                const Cell &stCell = fluid_->nearestItem(2 * cell.centroid() - nb.cell().centroid());
                LineSegment2D ln = intersectionLine(LineSegment2D(cell.centroid(), nb.cell().centroid()));
                Vector2D eta = nb.rCellVec().unitVec();

                Scalar eta1 = dot(stCell.centroid(), eta);
                Scalar eta2 = dot(cell.centroid(), eta);
                Scalar eta3 = dot(nb.cell().centroid(), eta);

                Matrix A = Matrix(3, 3, {
                        eta1 * eta1, eta1, 1.,
                        eta2 * eta2, eta2, 1.,
                        eta3 * eta3, eta3, 1.
                });

                Matrix c = solve(A, Matrix(3, 1, {
                        p(stCell),
                        p(cell),
                        p(nb.cell())
                }));

                Scalar etaB = dot(ln.ptB(), eta);

                bPts.push_back(ln.ptB());
                bP.push_back(c(0, 0) * etaB * etaB + c(1, 0) * etaB + c(2, 0) + rho * dot(g, ln.ptB()));
            }

    bPts = grid_.comm().allGatherv(bPts);
    bP = grid_.comm().allGatherv(bP);

    std::vector<ScalarPoint> pPoints;
    std::transform(bPts.begin(), bPts.end(), bP.begin(), std::back_inserter(pPoints), [](const Point2D &pt, Scalar p) {
        return std::make_pair(pt, p);
    });

    std::sort(pPoints.begin(), pPoints.end(), [this](const ScalarPoint &ptA, const ScalarPoint &ptB) {
        Vector2D rVecA = ptA.first - shapePtr_->centroid();
        Vector2D rVecB = ptB.first - shapePtr_->centroid();
        Scalar thetaA = std::atan2(rVecA.y, rVecA.x);
        Scalar thetaB = std::atan2(rVecB.y, rVecB.x);
        return (thetaA < 0. ? thetaA + 2 * M_PI : thetaA) < (thetaB < 0. ? thetaB + 2 * M_PI : thetaB);
    });

    force_ = Vector2D(0., 0.);

    for (int i = 0, end = pPoints.size(); i != end; ++i)
    {
        auto ptA = pPoints[i];
        auto ptB = pPoints[(i + 1) % end];
        Vector2D sf = (ptB.first - ptA.first).normalVec();
        force_ += -0.5 * (ptA.second + ptB.second) * sf;
    }

    //- weight
    force_ += this->rho * shapePtr_->area() * g;
}

void StepImmersedBoundaryObject::computeForce(const ScalarFiniteVolumeField &rho,
                                              const ScalarFiniteVolumeField &mu,
                                              const VectorFiniteVolumeField &u,
                                              const ScalarFiniteVolumeField &p,
                                              const Vector2D &g)
{
    typedef std::pair<Point2D, Scalar> ScalarPoint;
    typedef std::pair<Point2D, Vector2D> VectorPoint;

    std::vector<Point2D> bPts;
    std::vector<Scalar> bP;

    for (const Cell &cell: ibCells_)
        for (const InteriorLink &nb: cell.neighbours())
            if (isInIb(nb.cell()))
            {
                const Cell &stCell = fluid_->nearestItem(2 * cell.centroid() - nb.cell().centroid());
                LineSegment2D ln = intersectionLine(LineSegment2D(cell.centroid(), nb.cell().centroid()));
                Vector2D eta = nb.rCellVec().unitVec();

                Scalar eta1 = dot(stCell.centroid(), eta);
                Scalar eta2 = dot(cell.centroid(), eta);
                Scalar eta3 = dot(nb.cell().centroid(), eta);

                Matrix A = Matrix(3, 3, {
                        eta1 * eta1, eta1, 1.,
                        eta2 * eta2, eta2, 1.,
                        eta3 * eta3, eta3, 1.
                });

                Matrix c = solve(A, Matrix(3, 1, {
                        p(stCell),
                        p(cell),
                        p(nb.cell())
                }));

                Scalar etaB = dot(ln.ptB(), eta);
                Scalar pB = c(0, 0) * etaB * etaB + c(1, 0) * etaB + c(2, 0);

                A = Matrix(2, 2, {
                        eta1, 1.,
                        eta2, 1.
                });

                c = solve(A, Matrix(3, 1, {
                        rho(stCell),
                        rho(cell)
                }));

                Scalar rhoB = c(0, 0) * etaB + c(1, 0);

                bPts.push_back(ln.ptB());
                bP.push_back(pB + rhoB * dot(g, ln.ptB()));
            }

    bPts = grid_.comm().allGatherv(bPts);
    bP = grid_.comm().allGatherv(bP);

    std::vector<ScalarPoint> pPoints;
    std::transform(bPts.begin(), bPts.end(), bP.begin(), std::back_inserter(pPoints), [](const Point2D &pt, Scalar p) {
        return std::make_pair(pt, p);
    });

    std::sort(pPoints.begin(), pPoints.end(), [this](const ScalarPoint &ptA, const ScalarPoint &ptB) {
        Vector2D rVecA = ptA.first - shapePtr_->centroid();
        Vector2D rVecB = ptB.first - shapePtr_->centroid();
        Scalar thetaA = std::atan2(rVecA.y, rVecA.x);
        Scalar thetaB = std::atan2(rVecB.y, rVecB.x);
        return (thetaA < 0. ? thetaA + 2 * M_PI : thetaA) < (thetaB < 0. ? thetaB + 2 * M_PI : thetaB);
    });

    force_ = Vector2D(0., 0.);

    for (int i = 0, end = pPoints.size(); i != end; ++i)
    {
        auto ptA = pPoints[i];
        auto ptB = pPoints[(i + 1) % end];
        Vector2D sf = (ptB.first - ptA.first).normalVec();
        force_ += -0.5 * (ptA.second + ptB.second) * sf;
    }

    //- weight
    force_ += this->rho * shapePtr_->area() * g;
}
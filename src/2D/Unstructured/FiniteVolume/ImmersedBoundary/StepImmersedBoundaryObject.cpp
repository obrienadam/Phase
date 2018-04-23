#include <math.h>

#include "Math/Matrix.h"

#include "StepImmersedBoundaryObject.h"

void StepImmersedBoundaryObject::updateCells()
{
    solverCells_->add(cells_);
    auto cells = solverCells_->itemsWithin(*shape_);
    cells_.add(cells.begin(), cells.end());

    solidCells_.clear();
    solidCells_.add(cells_);

    ibCells_.clear();
    for (const Cell &cell: solidCells_)
        for (const InteriorLink &nb: cell.neighbours())
            if (!isInIb(nb.cell().centroid()))
                ibCells_.add(nb.cell()); //- Note: these should still be fluid cells
}

FiniteVolumeEquation<Scalar> StepImmersedBoundaryObject::bcs(ScalarFiniteVolumeField &field) const
{
    FiniteVolumeEquation<Scalar> eqn(field);

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
        case NORMAL_GRADIENT:
            for(const Cell& cell: cells_)
            {
                eqn.add(cell, cell, 1.);
                eqn.addSource(cell, -1.);
            }
    }

    return eqn;
}

FiniteVolumeEquation<Vector2D> StepImmersedBoundaryObject::bcs(VectorFiniteVolumeField &field) const
{
    FiniteVolumeEquation<Vector2D> eqn(field);

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
        case NORMAL_GRADIENT:
            for (const Cell &cell: cells_)
            {
                eqn.add(cell, cell, 1.);
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

    std::vector<Point2D> bPts, tauPtsX, tauPtsY;
    std::vector<Scalar> bP, tauX, tauY;

    for (const Cell &cell: ibCells_)
        for (const InteriorLink &nb: cell.neighbours())
            if (isInIb(nb.cell()))
            {
                const Cell &stCell = grid_->globalCells().nearestItem(
                        2. * cell.centroid() - nb.cell().centroid());

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

                if (std::abs(nb.rCellVec().x) > std::abs(nb.rCellVec().y))
                {
                    c = solve(A, Matrix(3, 1, {
                            u(stCell).y,
                            u(cell).y,
                            u(nb.cell()).y
                    }));

                    tauPtsY.push_back(ln.ptB());
                    tauY.push_back(-mu * (2. * c(0, 0) * etaB + c(1, 0)));
                }
                else
                {
                    c = solve(A, Matrix(3, 1, {
                            u(stCell).x,
                            u(cell).x,
                            u(nb.cell()).x
                    }));

                    tauPtsX.push_back(ln.ptB());
                    tauX.push_back(-mu * (2. * c(0, 0) * etaB + c(1, 0)));
                }
            }

    bPts = grid_->comm().allGatherv(bPts);
    bP = grid_->comm().allGatherv(bP);

    tauPtsX = grid_->comm().allGatherv(tauPtsX);
    tauX = grid_->comm().allGatherv(tauX);

    tauPtsY = grid_->comm().allGatherv(tauPtsY);
    tauY = grid_->comm().allGatherv(tauY);

    std::vector<ScalarPoint> pPoints, tauXPoints, tauYPoints;

    std::transform(bPts.begin(), bPts.end(), bP.begin(), std::back_inserter(pPoints), [](const Point2D &pt, Scalar p) {
        return std::make_pair(pt, p);
    });

    std::transform(tauPtsX.begin(), tauPtsX.end(), tauX.begin(), std::back_inserter(tauXPoints),
                   [](const Point2D &pt, Scalar p) {
                       return std::make_pair(pt, p);
                   });

    std::transform(tauPtsY.begin(), tauPtsY.end(), tauY.begin(), std::back_inserter(tauYPoints),
                   [](const Point2D &pt, Scalar p) {
                       return std::make_pair(pt, p);
                   });

    std::sort(pPoints.begin(), pPoints.end(), [this](const ScalarPoint &ptA, const ScalarPoint &ptB) {
        Vector2D rVecA = ptA.first - shape_->centroid();
        Vector2D rVecB = ptB.first - shape_->centroid();
        Scalar thetaA = std::atan2(rVecA.y, rVecA.x);
        Scalar thetaB = std::atan2(rVecB.y, rVecB.x);
        return (thetaA < 0. ? thetaA + 2. * M_PI : thetaA) < (thetaB < 0. ? thetaB + 2. * M_PI : thetaB);
    });

    std::sort(tauXPoints.end(), tauXPoints.end(), [this](const ScalarPoint &ptA, const ScalarPoint &ptB) {
        Vector2D rVecA = ptA.first - shape_->centroid();
        Vector2D rVecB = ptB.first - shape_->centroid();
        Scalar thetaA = std::atan2(rVecA.y, rVecA.x);
        Scalar thetaB = std::atan2(rVecB.y, rVecB.x);
        return (thetaA < 0. ? thetaA + 2. * M_PI : thetaA) < (thetaB < 0. ? thetaB + 2. * M_PI : thetaB);
    });

    std::sort(tauYPoints.end(), tauYPoints.end(), [this](const ScalarPoint &ptA, const ScalarPoint &ptB) {
        Vector2D rVecA = ptA.first - shape_->centroid();
        Vector2D rVecB = ptB.first - shape_->centroid();
        Scalar thetaA = std::atan2(rVecA.y, rVecA.x);
        Scalar thetaB = std::atan2(rVecB.y, rVecB.x);
        return (thetaA < 0. ? thetaA + 2. * M_PI : thetaA) < (thetaB < 0. ? thetaB + 2. * M_PI : thetaB);
    });

    force_ = Vector2D(0., 0.);
    torque_ = 0.;

    for (int i = 0, end = pPoints.size(); i != end; ++i)
    {
        auto ptA = pPoints[i];
        auto ptB = pPoints[(i + 1) % end];
        Vector2D sf = (ptB.first - ptA.first).normalVec();
        Vector2D xc = (ptA.first + ptB.first) / 2.;
        Vector2D df = -0.5 * (ptA.second + ptB.second) * sf;

        force_ += df;
        torque_ += cross(df, xc - shape_->centroid());
    }

    Vector2D shear(0., 0.);

    for (int i = 0, end = tauXPoints.size(); i != end; ++i)
    {
        auto ptA = tauXPoints[i];
        auto ptB = tauXPoints[(i + 1) % end];
        Vector2D sf = (ptB.first - ptA.first).normalVec();
        Vector2D xc = (ptA.first + ptB.first) / 2.;
        Vector2D df = Vector2D(-0.5 * (ptA.second + ptB.second) * sf.y, 0.);

        shear += df;
        torque_ += cross(df, xc - shape_->centroid());
    }

    for (int i = 0, end = tauYPoints.size(); i != end; ++i)
    {
        auto ptA = tauYPoints[i];
        auto ptB = tauYPoints[(i + 1) % end];
        Vector2D sf = (ptB.first - ptA.first).normalVec();
        Vector2D xc = (ptA.first + ptB.first) / 2.;
        Vector2D df = Vector2D(0., -0.5 * (ptA.second + ptB.second) * sf.x);

        shear += df;
        torque_ += cross(df, xc - shape_->centroid());
    }

    //- add weight and sheare to net force
    force_ += this->rho * shape_->area() * g + shear;
}

void StepImmersedBoundaryObject::computeForce(const ScalarFiniteVolumeField &rho,
                                              const ScalarFiniteVolumeField &mu,
                                              const VectorFiniteVolumeField &u,
                                              const ScalarFiniteVolumeField &p,
                                              const Vector2D &g)
{
    typedef std::pair<Point2D, Scalar> ScalarPoint;

    std::vector<Point2D> bPts;
    std::vector<Scalar> bP;

    for (const Cell &cell: ibCells_)
        for (const InteriorLink &nb: cell.neighbours())
            if (isInIb(nb.cell()))
            {
                const Cell &stCell = solverCells_->nearestItem(2 * cell.centroid() - nb.cell().centroid());
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

    bPts = grid_->comm().allGatherv(bPts);
    bP = grid_->comm().allGatherv(bP);

    std::vector<ScalarPoint> pPoints;
    std::transform(bPts.begin(), bPts.end(), bP.begin(), std::back_inserter(pPoints), [](const Point2D &pt, Scalar p) {
        return std::make_pair(pt, p);
    });

    std::sort(pPoints.begin(), pPoints.end(), [this](const ScalarPoint &ptA, const ScalarPoint &ptB) {
        Vector2D rVecA = ptA.first - shape_->centroid();
        Vector2D rVecB = ptB.first - shape_->centroid();
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
    force_ += this->rho * shape_->area() * g;
}
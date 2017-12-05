#include "QuadraticImmersedBoundaryObject.h"
#include "StaticMatrix.h"

QuadraticImmersedBoundaryObject::QuadraticImmersedBoundaryObject(const std::string &name,
                                                                 Label id,
                                                                 FiniteVolumeGrid2D &grid)
        :
        ImmersedBoundaryObject(name, id, grid)
{

}

//- Update
void QuadraticImmersedBoundaryObject::updateCells()
{
    fluid_->add(cells_);
    cells_.addAll(fluid_->itemsWithin(*shapePtr_));

    auto isIbCell = [this](const Cell &cell) {
        for (const InteriorLink &nb: cell.neighbours())
            if (!isInIb(nb.cell()))
                return true;

        for (const CellLink &dg: cell.diagonals())
            if (!isInIb(dg.cell()))
                return true;

        return false;
    };

    solidCells_.clear();
    ibCells_.clear();

    for (const Cell &cell: cells_)
    {
        if (isIbCell(cell))
            ibCells_.add(cell);
        else
            solidCells_.add(cell);
    }

    forcingCells_.clear();
    for (const Cell &cell: ibCells_)
        for (const InteriorLink &nb: cell.neighbours())
            if (!isInIb(nb.cell()))
                forcingCells_.add(nb.cell());
}

//- Boundary conditions
Equation <Scalar> QuadraticImmersedBoundaryObject::bcs(ScalarFiniteVolumeField &field) const
{
    Equation <Scalar> eqn(field);

    auto bType = boundaryType(field.name());
    auto bRefValue = getBoundaryRefValue<Scalar>(field.name());

    switch (bType)
    {
        case FIXED:
            for (const Cell &cell: solidCells_)
            {
                eqn.add(cell, cell, 1.);
                eqn.addSource(cell, -bRefValue);
            }

            for (const Cell &cell: ibCells_)
            {
                eqn.add(cell, cell, 1.);
                eqn.addSource(cell, -bRefValue);
            }

            break;
        default:
            throw Exception("QuadraticImmersedBoundaryObject", "bcs", "only fixed boundaries are supported.");
    }

    return eqn;
}

Equation <Vector2D> QuadraticImmersedBoundaryObject::bcs(VectorFiniteVolumeField &field) const
{
    Equation <Vector2D> eqn(field);

    auto bType = boundaryType(field.name());
    auto bRefValue = getBoundaryRefValue<Scalar>(field.name());

    switch (bType)
    {
        case FIXED:
            for (const Cell &cell: solidCells_)
            {
                eqn.add(cell, cell, 1.);
                eqn.addSource(cell, -bRefValue);
            }

            for (const Cell &cell: ibCells_)
            {
                eqn.add(cell, cell, 1.);
                eqn.addSource(cell, -bRefValue);
            }

            break;
        default:
            throw Exception("QuadraticImmersedBoundaryObject", "bcs", "only fixed boundaries are supported.");
    }

    return eqn;
}

Equation <Vector2D> QuadraticImmersedBoundaryObject::velocityBcs(VectorFiniteVolumeField &u) const
{
    Equation <Vector2D> eqn(u);

    auto bType = boundaryType(u.name());

    switch (bType)
    {
        case FIXED:
            for (const Cell &cell: solidCells_)
            {
                eqn.add(cell, cell, 1.);
                eqn.addSource(cell, -velocity(cell.centroid()));
            }

            for (const Cell &cell: ibCells_)
            {
                eqn.add(cell, cell, 1.);
                eqn.addSource(cell, -velocity(cell.centroid()));
            }

            break;
        default:
            throw Exception("QuadraticImmersedBoundaryObject", "velocityBcs", "only fixed boundaries are supported.");
    }

    return eqn;
}

void QuadraticImmersedBoundaryObject::computeForce(Scalar rho,
                                                   Scalar mu,
                                                   const VectorFiniteVolumeField &u,
                                                   const ScalarFiniteVolumeField &p,
                                                   const Vector2D &g)
{
    typedef std::pair <Point2D, Scalar> ScalarPoint;
    typedef std::pair <Point2D, Vector2D> VectorPoint;

    std::vector <Point2D> bPts, tauPtsX, tauPtsY;
    std::vector <Scalar> bP, tauX, tauY;
    bP.reserve(2 * ibCells_.size());

    for (const Cell &cell: ibCells_)
        for (const InteriorLink &nb: cell.neighbours())
            if (!isInIb(nb.cell()))
            {
                const Cell &cell0 = cell;
                const Cell &cell1 = nb.cell();
                const Cell &cell2 = grid_.globalActiveCells().nearestItem(2. * cell1.centroid() - cell0.centroid());

                if(cell1.id() == cell2.id())
                    continue;

                LineSegment2D ln = intersectionLine(LineSegment2D(cell0.centroid(), cell1.centroid()));
                Vector2D eta = (cell0.centroid() - cell1.centroid()).unitVec();

                Scalar eta0 = dot(cell0.centroid(), eta);
                Scalar eta1 = dot(cell1.centroid(), eta);
                Scalar eta2 = dot(cell2.centroid(), eta);



//                auto A = StaticMatrix<2, 2>(
//                        {
//                                eta1, 1.,
//                                eta2, 1.
//                        });
//
//                auto c = solve(A, StaticMatrix<2, 1>(
//                        {
//                                p(cell1),
//                                p(cell2)
//                        }));

                Scalar etaB = dot(ln.ptB(), eta);

                bPts.push_back(ln.ptB());
                bP.push_back((p(cell1) - p(cell2)) / (eta1 - eta2) * (etaB - eta2) + p(cell2) + rho * dot(ln.ptB(), g));

//                bP.push_back(c(0, 0) * etaB + c(1, 0) + rho * dot(g, ln.ptB()));


//                if (std::abs(nb.rCellVec().x) > std::abs(nb.rCellVec().y))
//                {
//                    c = solve(A, StaticMatrix<2, 1>(
//                            {
//                                    u(cell1).y,
//                                    u(cell2).y
//                            }));
//
//                    tauPtsY.push_back(ln.ptB());
//                    tauY.push_back(-mu * c(0, 0));
//                }
//                else
//                {
//                    c = solve(A, StaticMatrix<2, 1>(
//                            {
//                                    u(cell1).x,
//                                    u(cell2).x
//                            }));
//
//                    tauPtsX.push_back(ln.ptB());
//                    tauX.push_back(-mu * c(0, 0));
//                }
            }

    bPts = grid_.comm().allGatherv(bPts);
    bP = grid_.comm().allGatherv(bP);

    tauPtsX = grid_.comm().allGatherv(tauPtsX);
    tauX = grid_.comm().allGatherv(tauX);

    tauPtsY = grid_.comm().allGatherv(tauPtsY);
    tauY = grid_.comm().allGatherv(tauY);

    std::vector <ScalarPoint> pPoints, tauXPoints, tauYPoints;

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
        Vector2D rVecA = ptA.first - shapePtr_->centroid();
        Vector2D rVecB = ptB.first - shapePtr_->centroid();
        return rVecA.angle() < rVecB.angle();
    });

    std::sort(tauXPoints.end(), tauXPoints.end(), [this](const ScalarPoint &ptA, const ScalarPoint &ptB) {
        Vector2D rVecA = ptA.first - shapePtr_->centroid();
        Vector2D rVecB = ptB.first - shapePtr_->centroid();
        return rVecA.angle() < rVecB.angle();
    });

    std::sort(tauYPoints.end(), tauYPoints.end(), [this](const ScalarPoint &ptA, const ScalarPoint &ptB) {
        Vector2D rVecA = ptA.first - shapePtr_->centroid();
        Vector2D rVecB = ptB.first - shapePtr_->centroid();
        return rVecA.angle() < rVecB.angle();
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
        torque_ += cross(df, xc - shapePtr_->centroid());
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
        torque_ += cross(df, xc - shapePtr_->centroid());
    }

    for (int i = 0, end = tauYPoints.size(); i != end; ++i)
    {
        auto ptA = tauYPoints[i];
        auto ptB = tauYPoints[(i + 1) % end];
        Vector2D sf = (ptB.first - ptA.first).normalVec();
        Vector2D xc = (ptA.first + ptB.first) / 2.;
        Vector2D df = Vector2D(0., -0.5 * (ptA.second + ptB.second) * sf.x);

        shear += df;
        torque_ += cross(df, xc - shapePtr_->centroid());
    }

    //- add weight and sheare to net force
    force_ += this->rho * shapePtr_->area() * g;// + shear;
}

void QuadraticImmersedBoundaryObject::computeForce(const ScalarFiniteVolumeField &rho,
                                                   const ScalarFiniteVolumeField &mu,
                                                   const VectorFiniteVolumeField &u,
                                                   const ScalarFiniteVolumeField &p,
                                                   const Vector2D &g)
{
    typedef std::pair <Point2D, Scalar> ScalarPoint;
    typedef std::pair <Point2D, Vector2D> VectorPoint;

    std::vector <Point2D> bPts, tauPtsX, tauPtsY;
    std::vector <Scalar> bP, tauX, tauY;
    bP.reserve(2 * ibCells_.size());

    for (const Cell &cell: ibCells_)
        for (const InteriorLink &nb: cell.neighbours())
            if (!isInIb(nb.cell()))
            {
                const Cell &cell0 = cell;
                const Cell &cell1 = nb.cell();
                const Cell &cell2 = grid_.globalActiveCells().nearestItem(2. * cell1.centroid() - cell0.centroid());

                LineSegment2D ln = intersectionLine(LineSegment2D(cell1.centroid(), cell0.centroid()));
                Vector2D eta = (cell0.centroid() - cell1.centroid()).unitVec();

                Scalar eta1 = dot(cell1.centroid(), eta);
                Scalar eta2 = dot(cell2.centroid(), eta);
                Scalar etaB = dot(ln.ptB(), eta);

                Scalar pB = (p(cell1) - p(cell2)) / (eta1 - eta2) * (etaB - eta2) + p(cell2);
                Scalar rhoB = (rho(cell1) - rho(cell2)) / (eta1 - eta2) * (etaB - eta2) + rho(cell2);

                bPts.push_back(ln.ptB());
                bP.push_back(/*pB +*/ rhoB * dot(g, ln.ptB()));

//                if (std::abs(nb.rCellVec().x) > std::abs(nb.rCellVec().y))
//                {
//                    c = solve(A, StaticMatrix<2, 1>(
//                            {
//                                    u(cell1).y,
//                                    u(cell2).y
//                            }));
//
//                    Scalar dUdX = c(0, 0);
//
//                    c = solve(A, StaticMatrix<2, 1>(
//                            {
//                                    mu(cell1),
//                                    mu(cell2)
//                            }));
//
//                    Scalar muB = c(0, 0) * etaB + c(1, 0);
//
//                    tauPtsY.push_back(ln.ptB());
//                    tauY.push_back(-muB * dUdX);
//                }
//                else
//                {
//                    c = solve(A, StaticMatrix<2, 1>(
//                            {
//                                    u(cell1).x,
//                                    u(cell2).x
//                            }));
//
//                    Scalar dUdY = c(0, 0);
//
//                    c = solve(A, StaticMatrix<2, 1>(
//                            {
//                                    mu(cell1),
//                                    mu(cell2)
//                            }));
//
//                    Scalar muB = c(0, 0) * etaB + c(1, 0);
//                    muB = std::max(std::min(muB, 8.9e-4), 1.81e-5);
//
//                    tauPtsX.push_back(ln.ptB());
//                    tauX.push_back(-muB * dUdY);
//                }
            }

    bPts = grid_.comm().allGatherv(bPts);
    bP = grid_.comm().allGatherv(bP);

    tauPtsX = grid_.comm().allGatherv(tauPtsX);
    tauX = grid_.comm().allGatherv(tauX);

    tauPtsY = grid_.comm().allGatherv(tauPtsY);
    tauY = grid_.comm().allGatherv(tauY);

    std::vector <ScalarPoint> pPoints, tauXPoints, tauYPoints;
    pPoints.reserve(bP.size());

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
        Vector2D rVecA = ptA.first - shapePtr_->centroid();
        Vector2D rVecB = ptB.first - shapePtr_->centroid();
        return rVecA.angle() < rVecB.angle();
    });

    std::sort(tauXPoints.end(), tauXPoints.end(), [this](const ScalarPoint &ptA, const ScalarPoint &ptB) {
        Vector2D rVecA = ptA.first - shapePtr_->centroid();
        Vector2D rVecB = ptB.first - shapePtr_->centroid();
        return rVecA.angle() < rVecB.angle();
    });

    std::sort(tauYPoints.end(), tauYPoints.end(), [this](const ScalarPoint &ptA, const ScalarPoint &ptB) {
        Vector2D rVecA = ptA.first - shapePtr_->centroid();
        Vector2D rVecB = ptB.first - shapePtr_->centroid();
        return rVecA.angle() < rVecB.angle();
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
        torque_ += cross(df, xc - shapePtr_->centroid());
        Vector2D r = ptA.first - shapePtr_->centroid();
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
        torque_ += cross(df, xc - shapePtr_->centroid());
    }

    for (int i = 0, end = tauYPoints.size(); i != end; ++i)
    {
        auto ptA = tauYPoints[i];
        auto ptB = tauYPoints[(i + 1) % end];
        Vector2D sf = (ptB.first - ptA.first).normalVec();
        Vector2D xc = (ptA.first + ptB.first) / 2.;
        Vector2D df = Vector2D(0., -0.5 * (ptA.second + ptB.second) * sf.x);

        shear += df;
        torque_ += cross(df, xc - shapePtr_->centroid());
    }

    //- add weight and sheare to net force
    force_ += this->rho * shapePtr_->area() * g; // + shear;
}
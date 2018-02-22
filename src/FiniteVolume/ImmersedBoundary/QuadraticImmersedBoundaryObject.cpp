#include "QuadraticImmersedBoundaryObject.h"
#include "BilinearInterpolator.h"

QuadraticImmersedBoundaryObject::QuadraticImmersedBoundaryObject(const std::string &name,
                                                                 Label id,
                                                                 const ImmersedBoundary &ib,
                                                                 const std::shared_ptr<FiniteVolumeGrid2D> &grid)
        :
        ImmersedBoundaryObject(name, id, ib, grid)
{

}

//- Clear

void QuadraticImmersedBoundaryObject::clear()
{
    ImmersedBoundaryObject::clear();
    stencils_.clear();
    forcingCells_.clear();
}

//- Update
void QuadraticImmersedBoundaryObject::updateCells()
{
    //std::cout << "Add cells back to fluid for \"" << name_ << "\"...\n";
    clear();

    //std::cout << "Find cells inside \"" << name_ << "\"...\n";
    auto items = fluid_->itemsWithin(*shape_);

    //std::cout << "Add items to cell zone for \"" << name_ << "\"...\n";
    cells_.add(items.begin(), items.end());

    //std::cout << "Constructing cell zones for \"" << name_ << "\"...\n";
    auto isIbCell = [this](const Cell &cell)
    {
        for (const CellLink &nb: cell.neighbours())
            if (!isInIb(nb.cell()))
                return true;

        for (const CellLink &dg: cell.diagonals())
            if (!isInIb(dg.cell()))
                return true;

        return false;
    };

    for (const Cell &cell: cells_)
    {
        if (isIbCell(cell))
            ibCells_.add(cell);
        else
            solidCells_.add(cell);
    }

    for (const Cell &cell: ibCells_)
        for (const InteriorLink &nb: cell.neighbours())
            if (fluid_->isInGroup(nb.cell()))
                forcingCells_.add(nb.cell());

    for (const Cell &cell: forcingCells_)
        for (const InteriorLink &nb: cell.neighbours())
        {
            auto ibObj = ib_->ibObj(nb.cell().centroid());

            if (ibObj)
                stencils_.push_back(QuadraticIbmStencil(nb, *ib_));
        }
}

//- Boundary conditions
Equation<Scalar> QuadraticImmersedBoundaryObject::bcs(ScalarFiniteVolumeField &field) const
{
    Equation<Scalar> eqn(field);

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

Equation<Vector2D> QuadraticImmersedBoundaryObject::bcs(VectorFiniteVolumeField &field) const
{
    Equation<Vector2D> eqn(field);

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

Equation<Vector2D> QuadraticImmersedBoundaryObject::velocityBcs(VectorFiniteVolumeField &u) const
{
    Equation<Vector2D> eqn(u);

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
    std::vector<Point2D> points;
    std::vector<Scalar> pressures;
    std::vector<Scalar> shears;

    points.reserve(ibCells_.size());
    pressures.reserve(ibCells_.size());
    shears.reserve(ibCells_.size());

    auto bi = BilinearInterpolator(grid_);
    for (const Cell &cell: ibCells_)
    {
        Point2D pt = nearestIntersect(cell.centroid());
        Vector2D wn = nearestEdgeNormal(pt).unitVec();
        bi.setPoint(pt);

        if (bi.isValid())
        {
            points.push_back(pt);
            pressures.push_back(bi(p) + rho * dot(pt, g));
            shears.push_back(mu * dot(dot(bi.grad(u), wn), wn.tangentVec()));
        }
    }

    points = grid_->comm().gatherv(grid_->comm().mainProcNo(), points);
    pressures = grid_->comm().gatherv(grid_->comm().mainProcNo(), pressures);
    shears = grid_->comm().gatherv(grid_->comm().mainProcNo(), shears);

    if (grid_->comm().isMainProc())
    {
        force_ = Vector2D(0., 0.);

        std::vector<std::tuple<Point2D, Scalar, Scalar>> stresses(points.size());
        for (int i = 0; i < points.size(); ++i)
            stresses[i] = std::make_tuple(points[i], pressures[i], shears[i]);

        std::sort(stresses.begin(), stresses.end(),
                  [this](const std::tuple<Point2D, Scalar, Scalar> &a, std::tuple<Point2D, Scalar, Scalar> &b)
                  {
                      return (std::get<0>(a) - shape_->centroid()).angle() <
                             (std::get<0>(b) - shape_->centroid()).angle();
                  });

        for (int i = 0; i < stresses.size(); ++i)
        {
            const auto &a = stresses[i];
            const auto &b = stresses[(i + 1) % stresses.size()];

            const Point2D &ptA = std::get<0>(a);
            const Point2D &ptB = std::get<0>(b);
            Scalar prA = std::get<1>(a);
            Scalar prB = std::get<1>(b);
            Scalar shA = std::get<2>(a);
            Scalar shB = std::get<2>(b);

            force_ += -(prA + prB) / 2. * (ptB - ptA).normalVec() + (shA + shB) / 2. * (ptB - ptA);
        }
    }

    force_ = grid_->comm().broadcast(grid_->comm().mainProcNo(), force_) + shape_->area() * this->rho * g;

    auto maxF = grid_->comm().max(force_.mag());

    if (grid_->comm().isMainProc())
        std::cout << "MAX PARTICLE FORCE = " << maxF << std::endl;
}

void QuadraticImmersedBoundaryObject::computeForce(const ScalarFiniteVolumeField &rho,
                                                   const ScalarFiniteVolumeField &mu,
                                                   const VectorFiniteVolumeField &u,
                                                   const ScalarFiniteVolumeField &p,
                                                   const Vector2D &g)
{
    std::vector<Point2D> points;
    std::vector<Scalar> pressures;
    std::vector<Scalar> shears;
    points.reserve(ibCells_.size());
    pressures.reserve(ibCells_.size());
    shears.reserve(ibCells_.size());

    auto bi = BilinearInterpolator(grid_);
    for (const Cell &cell: ibCells_)
    {
        Point2D pt = nearestIntersect(cell.centroid());
        Vector2D wn = nearestEdgeNormal(pt).unitVec();
        bi.setPoint(pt);

        if (bi.isValid())
        {
            points.push_back(pt);
            pressures.push_back(bi(p) + bi(rho) * (dot(pt, g) + 0.5 * std::pow(dot(bi(u), wn), 2)));
            shears.push_back(bi(mu) * dot(dot(bi.grad(u), wn), wn.tangentVec()));
        }
    }

    points = grid_->comm().gatherv(grid_->comm().mainProcNo(), points);
    pressures = grid_->comm().gatherv(grid_->comm().mainProcNo(), pressures);
    shears = grid_->comm().gatherv(grid_->comm().mainProcNo(), shears);

    if (grid_->comm().isMainProc())
    {
        force_ = Vector2D(0., 0.);

        std::vector<std::tuple<Point2D, Scalar, Scalar>> stresses(points.size());

        for (int i = 0; i < points.size(); ++i)
            stresses[i] = std::make_tuple(points[i], pressures[i], shears[i]);

        std::sort(stresses.begin(), stresses.end(),
                  [this](const std::tuple<Point2D, Scalar, Scalar> &a,
                         const std::tuple<Point2D, Scalar, Scalar> &b)
                  {
                      return (std::get<0>(a) - shape_->centroid()).angle() <
                             (std::get<0>(b) - shape_->centroid()).angle();
                  });

        for (int i = 0; i < stresses.size(); ++i)
        {
            const auto &a = stresses[i];
            const auto &b = stresses[(i + 1) % stresses.size()];

            const Point2D &ptA = std::get<0>(a);
            const Point2D &ptB = std::get<0>(b);
            Scalar prA = std::get<1>(a);
            Scalar prB = std::get<1>(b);
            Scalar shA = std::get<2>(a);
            Scalar shB = std::get<2>(b);

            force_ += -(prA + prB) / 2. * (ptB - ptA).normalVec() + (shA + shB) / 2. * (ptB - ptA);
        }
    }

    force_ = grid_->comm().broadcast(grid_->comm().mainProcNo(), force_) + shape_->area() * this->rho * g;
}

void QuadraticImmersedBoundaryObject::computeForce(Scalar rho1,
                                                   Scalar rho2,
                                                   Scalar mu1,
                                                   Scalar mu2,
                                                   const ScalarFiniteVolumeField &gamma,
                                                   const VectorFiniteVolumeField &u, const ScalarFiniteVolumeField &p,
                                                   const Vector2D &g)
{

}
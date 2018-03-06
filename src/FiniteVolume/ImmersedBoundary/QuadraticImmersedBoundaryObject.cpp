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
    clear();
    auto items = fluid_->itemsCoveredBy(*shape_);
    solidCells_.add(items.begin(), items.end());
    cells_.add(solidCells_);
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
        default:throw Exception("QuadraticImmersedBoundaryObject", "bcs", "only fixed boundaries are supported.");
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
        default:throw Exception("QuadraticImmersedBoundaryObject", "bcs", "only fixed boundaries are supported.");
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
    std::vector<std::tuple<Point2D, Scalar, Scalar>> stresses;
    stresses.reserve(ibCells_.size());

    auto bi = BilinearInterpolator(grid_);
    for (const Cell &cell: ibCells_)
    {
        Point2D pt = nearestIntersect(cell.centroid());
        Vector2D wn = nearestEdgeNormal(pt).unitVec();
        bi.setPoint(pt);

        if (bi.isValid())
        {
            stresses.push_back(
                    std::make_tuple(
                            pt,
                            bi(p) + rho * (dot(pt, g) + 0.5 * std::pow(dot(bi(u), wn), 2)),
                            mu * dot(dot(bi.grad(u), wn), wn.tangentVec())
                    )
            );
        }
    }

    stresses = grid_->comm().gatherv(grid_->comm().mainProcNo(), stresses);

    if (grid_->comm().isMainProc())
    {
        std::sort(stresses.begin(), stresses.end(),
                  [this](const std::tuple<Point2D, Scalar, Scalar> &a,
                         const std::tuple<Point2D, Scalar, Scalar> &b)
                  {
                      return (std::get<0>(a) - shape_->centroid()).angle() <
                             (std::get<0>(b) - shape_->centroid()).angle();
                  });

        force_ = Vector2D(0., 0.);

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

void QuadraticImmersedBoundaryObject::computeForce(const ScalarFiniteVolumeField &rho,
                                                   const ScalarFiniteVolumeField &mu,
                                                   const VectorFiniteVolumeField &u,
                                                   const ScalarFiniteVolumeField &p,
                                                   const Vector2D &g)
{
    std::vector<std::tuple<Point2D, Scalar, Scalar>> stresses;
    stresses.reserve(ibCells_.size());

    auto bi = BilinearInterpolator(grid_);
    for (const Cell &cell: ibCells_)
    {
        Point2D pt = nearestIntersect(cell.centroid());
        Vector2D wn = nearestEdgeNormal(pt).unitVec();
        bi.setPoint(pt);

        if (bi.isValid())
        {
            stresses.push_back(
                    std::make_tuple(
                            pt,
                            bi(p) + bi(rho) * dot(pt, g),
                            bi(mu) * dot(dot(bi.grad(u), wn), wn.tangentVec())
                    )
            );
        }
    }

    stresses = grid_->comm().gatherv(grid_->comm().mainProcNo(), stresses);

    if (grid_->comm().isMainProc())
    {
        std::sort(stresses.begin(), stresses.end(),
                  [this](const std::tuple<Point2D, Scalar, Scalar> &a,
                         const std::tuple<Point2D, Scalar, Scalar> &b)
                  {
                      return (std::get<0>(a) - shape_->centroid()).angle() <
                             (std::get<0>(b) - shape_->centroid()).angle();
                  });

        force_ = Vector2D(0., 0.);

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

void QuadraticImmersedBoundaryObject::computeForce(const ScalarFiniteVolumeField &rho,
                                                   const ScalarFiniteVolumeField &mu,
                                                   const VectorFiniteVolumeField &u,
                                                   const ScalarFiniteVolumeField &p,
                                                   const ScalarFiniteVolumeField &gamma,
                                                   const SurfaceTensionForce &ft,
                                                   const Vector2D &g)
{
    computeForce(rho, mu, u, p, g);
    force_ += ft.computeCapillaryForce(gamma, *this);
}
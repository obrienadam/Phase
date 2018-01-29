#include <tuple>

#include "GhostCellImmersedBoundaryObject.h"

GhostCellImmersedBoundaryObject::GhostCellImmersedBoundaryObject(const std::string &name,
                                                                 Label id,
                                                                 const ImmersedBoundary &ib,
                                                                 const std::shared_ptr<FiniteVolumeGrid2D> &grid)
        :
        ImmersedBoundaryObject(name, id, ib, grid)
{

}

void GhostCellImmersedBoundaryObject::updateCells()
{
    clear();

    auto items = fluid_->itemsWithin(*shape_);
    cells_.add(items.begin(), items.end());

    auto isIbCell = [this](const Cell &cell) {
        if (!isInIb(cell.centroid()))
            return false;

        for (const InteriorLink &nb: cell.neighbours())
            if (!isInIb(nb.cell().centroid()))
                return true;

        for (const CellLink &dg: cell.diagonals())
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

    switch (bType)
    {
        case FIXED:
            for (const GhostCellStencil &st: stencils_)
            {
                eqn.add(st.cell(), st.dirichletCells(), st.dirichletCoeffs());
                eqn.addSource(st.cell(), -bRefValue);
            }

            for (const Cell &cell: solidCells_)
            {
                eqn.add(cell, cell, 1.);
                eqn.addSource(cell, -bRefValue);
            }

            break;
        case NORMAL_GRADIENT:
            for (const GhostCellStencil &st: stencils_)
            {
                eqn.add(st.cell(), st.neumannCells(), st.neumannCoeffs());
                eqn.addSource(st.cell(), -bRefValue);
            }

            for (const Cell &cell: solidCells_)
                eqn.add(cell, cell, 1.);
            break;

        default:
            throw Exception("GhostCellImmersedBoundaryObject", "bcs", "invalid boundary type.");
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
                eqn.add(st.cell(), st.dirichletCells(), st.dirichletCoeffs());
                eqn.addSource(st.cell(), -bRefValue);
                break;

            case ImmersedBoundaryObject::NORMAL_GRADIENT:
                eqn.add(st.cell(), st.neumannCells(), st.neumannCoeffs());
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

Equation<Vector2D> GhostCellImmersedBoundaryObject::velocityBcs(VectorFiniteVolumeField &u) const
{
    Equation<Vector2D> eqn(u);

    switch (boundaryType(u.name()))
    {
        case FIXED:
            for (const GhostCellStencil &st: stencils_)
            {
                eqn.add(st.cell(), st.dirichletCells(), st.dirichletCoeffs());
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

Equation<Scalar> GhostCellImmersedBoundaryObject::pressureBcs(Scalar rho, ScalarFiniteVolumeField &p) const
{
    Equation<Scalar> eqn(p);

    for (const GhostCellStencil &st: stencils_)
        eqn.add(st.cell(), st.neumannCells(), st.neumannCoeffs());

    if (motion_)
        for (const GhostCellStencil &st: stencils_)
        {
            Scalar dUdN = dot(acceleration(st.boundaryPoint()), nearestEdgeNormal(st.boundaryPoint()).unitVec());
            eqn.addSource(st.cell(), rho * dUdN);
        }

    for (const Cell &cell: solidCells())
    {
        eqn.add(cell, cell, 1.);
    }

    return eqn;
}

Equation<Scalar> GhostCellImmersedBoundaryObject::contactLineBcs(ScalarFiniteVolumeField &gamma, Scalar theta) const
{
    Equation<Scalar> eqn(gamma);
    Scalar minDistSqr = 100;
    int i = 0;

    for (const GhostCellStencil &st: stencils_)
    {
        Vector2D wn = -nearestEdgeNormal(st.boundaryPoint());

        Ray2D r1 = Ray2D(st.cell().centroid(), wn.rotate(M_PI_2 - theta));
        Ray2D r2 = Ray2D(st.cell().centroid(), wn.rotate(theta - M_PI_2));

        GhostCellStencil m1(st.cell(), shape().intersections(r1)[0], r1.r(), *grid_);
        GhostCellStencil m2(st.cell(), shape().intersections(r2)[0], r2.r(), *grid_);
        Scalar g1 = m1.bpValue(gamma);
        Scalar g2 = m2.bpValue(gamma);

        if (std::abs(g1 - g2) > 1e-8)
        {
            if (g2 < g1)
                std::swap(m1, m2);
        }
        else
        {
            Vector2D grad1 = m1.bpGrad(gamma);
            Vector2D grad2 = m2.bpGrad(gamma);

            if(dot(grad1, r1.r()) < 0 && theta > M_PI_2)
                std::swap(m1, m2);
        }



//        if(!(g1 < 0.5 && g2 > 0.5))
//        {
//            eqn.add(st.cell(), st.cells(), st.neumannCoeffs());
//            continue;
//        }

        if (theta > M_PI_2)
            eqn.add(m1.cell(), m1.neumannCells(), m1.neumannCoeffs());
        else
            eqn.add(m2.cell(), m2.neumannCells(), m2.neumannCoeffs());

//        Scalar distSqr = (st.cell().centroid() - nearestIntersect(st.cell().centroid())).magSqr();
//
//        if (distSqr < minDistSqr && st.cell().centroid().y > 0.85)
//        {
//            minDistSqr = distSqr;
//            std::cout << "NEW MIN CELL FOUND.\n"
//                      << st.cell().centroid() << "\n"
//                      << m1.boundaryPoint() << "\n"
//                      << m2.boundaryPoint() << "\n"
//                      << m1.bpGrad(gamma) << "\n" << m2.bpGrad(gamma) << "\n"
//                      << g1 << "\n" << g2 << "\n";
//        }
    }

    for (const Cell &cell: solidCells_)
        eqn.add(cell, cell, 1.);

    return eqn;
}

//- Force computation

void GhostCellImmersedBoundaryObject::computeForce(Scalar rho,
                                                   Scalar mu,
                                                   const VectorFiniteVolumeField &u,
                                                   const ScalarFiniteVolumeField &p,
                                                   const Vector2D &g)
{
    std::vector<Point2D> points;
    std::vector<Scalar> pressures;
    std::vector<Scalar> shears;

    for (const GhostCellStencil &st: stencils_)
    {
        points.push_back(st.boundaryPoint());
        pressures.push_back(st.bpValue(p) + rho * dot(st.boundaryPoint(), g));
        shears.push_back(mu * dot(dot(st.bpGrad(u), st.wallNormal()), st.wallNormal().tangentVec()));
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
                  [this](const std::tuple<Point2D, Scalar, Scalar> &a, std::tuple<Point2D, Scalar, Scalar> &b) {
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

    force_ = grid_->comm().broadcast(grid_->comm().mainProcNo(), force_);
}

void GhostCellImmersedBoundaryObject::computeForce(const ScalarFiniteVolumeField &rho,
                                                   const ScalarFiniteVolumeField &mu,
                                                   const VectorFiniteVolumeField &u,
                                                   const ScalarFiniteVolumeField &p,
                                                   const Vector2D &g)
{

}

//- Protected methods

void GhostCellImmersedBoundaryObject::constructStencils()
{
    stencils_.clear();
    for (const Cell &cell: ibCells_)
        stencils_.push_back(GhostCellStencil(cell, *this));
}
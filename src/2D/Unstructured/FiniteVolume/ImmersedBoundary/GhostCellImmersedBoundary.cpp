#include <tuple>

#include "GhostCellImmersedBoundary.h"
#include "GhostCellImmersedBoundaryFixedBcStencil.h"
#include "GhostCellImmersedBoundaryNormalGradientBcStencil.h"

GhostCellImmersedBoundary::GhostCellImmersedBoundary(const Input &input,
                                                     const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                                     const std::shared_ptr<CellGroup> &domainCells)
    :
      ImmersedBoundary(input, grid, domainCells)
{

}

void GhostCellImmersedBoundary::updateCells()
{
    for(auto &ibObj: ibObjs_)
    {
        domainCells_->add(ibObj->cells());
        ibObj->clear();
    }

    CellGroup ibCells, solidCells;

    for(auto &ibObj: ibObjs_)
    {
        solidCells = ibObj->cellsWithin(*domainCells());
        ibCells = ibObj->internalPerimeterCells(solidCells.begin(), solidCells.end(), true);

        solidCells.remove(ibCells);

        domainCells_->remove(ibCells);
        domainCells_->remove(solidCells);
    }
}

FiniteVolumeEquation<Scalar> GhostCellImmersedBoundary::bcs(ScalarFiniteVolumeField &phi) const
{
    FiniteVolumeEquation<Scalar> eqn(phi);

    for(const auto &ibObj: ibObjs_)
    {
        Scalar refVal = ibObj->bcRefValue<Scalar>(phi.name());

        switch(ibObj->bcType(phi.name()))
        {
        case ImmersedBoundaryObject::FIXED:
            for(const Cell& cell: ibObj->ibCells())
            {
                FixedBcStencil st = FixedBcStencil(cell, *ibObj, *grid_);
                eqn.add(st.cell(), st.cells(), st.coeffs());
                eqn.addSource(st.cell(), -refVal);
            }

            for(const Cell& cell: ibObj->solidCells())
            {
                eqn.add(cell, cell, 1.);
                eqn.addSource(cell, -refVal);
            }

            break;
        case ImmersedBoundaryObject::NORMAL_GRADIENT:
            for(const Cell& cell: ibObj->ibCells())
            {
                NormalGradientBcStencil st = NormalGradientBcStencil(cell, *ibObj, *grid_);
                eqn.add(st.cell(), st.cells(), st.coeffs());
                eqn.addSource(st.cell(), -refVal);
            }

            for(const Cell& cell: ibObj->solidCells())
            {
                eqn.add(cell, cell, 1.);
                eqn.setSource(cell, 0.);
            }

            break;
        default:
            throw Exception("GhostCellImmersedBoundary", "bcs", "invalid boundary condition type.");
        }
    }

    return eqn;
}

FiniteVolumeEquation<Vector2D> GhostCellImmersedBoundary::velocityBcs(VectorFiniteVolumeField &u) const
{
    FiniteVolumeEquation<Vector2D> eqn(u);

    for(const auto &ibObj: ibObjs_)
    {
        for(const Cell& cell: ibObj->ibCells())
        {
            FixedBcStencil st = FixedBcStencil(cell, *ibObj, *grid_);
            eqn.add(st.cell(), st.cells(), st.coeffs());
            eqn.addSource(st.cell(), -ibObj->velocity(st.bp()));
        }

        for(const Cell& cell: ibObj->solidCells())
        {
            eqn.add(cell, cell, 1.);
            eqn.addSource(cell, -ibObj->velocity(cell.centroid()));
        }
    }

    return eqn;
}

//- Force computation

void GhostCellImmersedBoundary::applyHydrodynamicForce(Scalar rho,
                                                       Scalar mu,
                                                       const VectorFiniteVolumeField &u,
                                                       const ScalarFiniteVolumeField &p,
                                                       const Vector2D &g)
{
    std::vector<std::tuple<Point2D, Scalar, Scalar>> stresses;

    for(auto &ibObj: ibObjs_)
    {
        for (const Cell& cell: ibObj->ibCells())
        {
            FixedBcStencil st = FixedBcStencil(cell, *ibObj, *grid_);

            stresses.push_back(
                        std::make_tuple(
                            st.bp(),
                            st.bpValue(p) + rho * dot(st.bp(), g),
                            mu * dot(dot(st.bpGrad(u), st.nw()), st.nw().tangentVec())
                            )
                        );
        }

        stresses = grid_->comm().gatherv(grid_->comm().mainProcNo(), stresses);

        auto force = Vector2D(0., 0.);

        if (grid_->comm().isMainProc())
        {
            std::sort(stresses.begin(), stresses.end(),
                      [ibObj](const std::tuple<Point2D, Scalar, Scalar> &a, std::tuple<Point2D, Scalar, Scalar> &b)
            {
                return (std::get<0>(a) - ibObj->shape().centroid()).angle() <
                        (std::get<0>(b) - ibObj->shape().centroid()).angle();
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

                force += -(prA + prB) / 2. * (ptB - ptA).normalVec() + (shA + shB) / 2. * (ptB - ptA);
            }
        }

        force = grid_->comm().broadcast(grid_->comm().mainProcNo(), force);
        ibObj->applyForce(force);

        std::cout << "Force = " << force << std::endl;
    }
}

void GhostCellImmersedBoundary::applyHydrodynamicForce(const ScalarFiniteVolumeField &rho,
                                                       const ScalarFiniteVolumeField &mu,
                                                       const VectorFiniteVolumeField &u,
                                                       const ScalarFiniteVolumeField &p,
                                                       const Vector2D &g)
{
    throw Exception("GhostCellImmersedBoundary", "applyHydrodynamicForce", "not implemented.");
}

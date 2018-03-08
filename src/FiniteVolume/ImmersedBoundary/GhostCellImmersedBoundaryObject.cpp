#include <tuple>

#include "GhostCellImmersedBoundaryObject.h"

GhostCellImmersedBoundaryObject::Stencil::Stencil(const Cell &cell,
                                                  const ImmersedBoundaryObject &ibObj,
                                                  const Vector2D &d)
        :
        cell_(std::cref(cell)),
        d_(d)
{
    bp_ = ibObj.shape().intersections(Ray2D(cell.centroid(), d_))[0];
    ip_ = 2. * bp_ - cell.centroid();
    nw_ = -ibObj.nearestEdgeNormal(bp_);

    ipCells_ = ibObj.grid()->findNearestNode(ip_).cells();
    bpCells_ = ibObj.grid()->findNearestNode(bp_).cells();

    for (int i = 0; i < 4; ++i)
    {
        const Point2D &xi = ipCells_[i].get().centroid();
        const Point2D &xb = bpCells_[i].get().centroid();

        Aip_.setRow(i, {xi.x * xi.y, xi.x, xi.y, 1.});
        Abp_.setRow(i, {xb.x * xb.y, xb.x, xb.y, 1.});
    }

    Aip_.invert();
    Abp_.invert();
}

GhostCellImmersedBoundaryObject::Stencil::Stencil(const Cell &cell,
                                                  const ImmersedBoundaryObject &ibObj)
        :
        Stencil(cell, ibObj, -ibObj.nearestEdgeNormal(cell.centroid()))
{

}

void GhostCellImmersedBoundaryObject::Stencil::initFixedBc()
{
    fixedBcCells_ = ipCells_;
    fixedBcCells_.push_back(std::cref(cell_));
    auto cd = StaticMatrix<1, 4>({ip_.x * ip_.y, ip_.x, ip_.y, 1.}) * Aip_ / 2.;
    fixedBcCoeffs_.assign(cd.begin(), cd.end());
    fixedBcCoeffs_.push_back(0.5);
}

void GhostCellImmersedBoundaryObject::Stencil::initZeroGradientBc()
{
    zeroGradBcCells_ = bpCells_;
    zeroGradBcCells_.push_back(std::cref(cell_));

    bool bpInStencil = false;
    for (const Cell &cell: bpCells_)
        if (cell.id() == cell.id())
        {
            bpInStencil = true;
            break;
        }

    if (bpInStencil)
    {
        auto cn = StaticMatrix<1, 4>({bp_.y * d_.x + bp_.x * d_.y, d_.x, d_.y, 0.}) * Abp_;
        zeroGradBcCoeffs_.assign(cn.begin(), cn.end());
        zeroGradBcCoeffs_.push_back(0.);
    }
    else
    {
        auto cn = StaticMatrix<1, 4>({ip_.x * ip_.y, ip_.x, ip_.y, 1.}) * Abp_;
        zeroGradBcCoeffs_.assign(cn.begin(), cn.end());
        zeroGradBcCoeffs_.push_back(-1.);
    }
}

Scalar GhostCellImmersedBoundaryObject::Stencil::bpValue(const ScalarFiniteVolumeField &phi) const
{
    return (StaticMatrix<1, 4>({bp_.x * bp_.y, bp_.x, bp_.y, 1.})
            * Abp_
            * StaticMatrix<4, 1>(
                    {
                            phi(bpCells_[0]),
                            phi(bpCells_[1]),
                            phi(bpCells_[2]),
                            phi(bpCells_[3])
                    }))(0, 0);
}

Vector2D GhostCellImmersedBoundaryObject::Stencil::bpValue(const VectorFiniteVolumeField &u) const
{
    auto ubp = StaticMatrix<1, 4>({bp_.x * bp_.y, bp_.x, bp_.y, 1.})
               * Abp_
               * StaticMatrix<4, 2>(
                       {
                               u(bpCells_[0]).x, u(bpCells_[0]).y,
                               u(bpCells_[1]).x, u(bpCells_[1]).y,
                               u(bpCells_[2]).x, u(bpCells_[2]).y,
                               u(bpCells_[3]).x, u(bpCells_[3]).y,
                       });

    return Vector2D(ubp(0, 0), ubp(0, 1));
}

Vector2D GhostCellImmersedBoundaryObject::Stencil::bpGrad(const ScalarFiniteVolumeField &phi) const
{
    auto x = StaticMatrix<2, 4>(
            {
                    bp_.y, 1., 0., 0.,
                    bp_.x, 0., 1., 0.
            })
             * Abp_
             * StaticMatrix<4, 1>(
                     {
                             phi(bpCells_[0]),
                             phi(bpCells_[1]),
                             phi(bpCells_[2]),
                             phi(bpCells_[3]),
                     });

    return Vector2D(x(0, 0), x(1, 0));
}

Tensor2D GhostCellImmersedBoundaryObject::Stencil::bpGrad(const VectorFiniteVolumeField &u) const
{
    auto x = StaticMatrix<2, 4>(
            {
                    bp_.y, 1., 0., 0.,
                    bp_.x, 0., 1., 0.
            })
             * Abp_
             * StaticMatrix<4, 2>(
                     {
                             u(bpCells_[0]).x, u(bpCells_[0]).y,
                             u(bpCells_[1]).x, u(bpCells_[1]).y,
                             u(bpCells_[2]).x, u(bpCells_[2]).y,
                             u(bpCells_[3]).x, u(bpCells_[3]).y
                     });

    return Tensor2D(x(0, 0), x(0, 1), x(1, 0), x(1, 1));
}

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

    auto isIbCell = [this](const Cell &cell)
    {
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

    fluid_->add(solidCells_);
    solidCells_.clear();
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
            for (const Stencil &st: stencils_)
            {
                eqn.add(st.cell(), st.fixedBcCells(), st.fixedBcCoeffs());
                eqn.addSource(st.cell(), -bRefValue);
            }

            for (const Cell &cell: solidCells_)
            {
                eqn.add(cell, cell, 1.);
                eqn.addSource(cell, -bRefValue);
            }

            break;
        case NORMAL_GRADIENT:
            for (const Stencil &st: stencils_)
                eqn.add(st.cell(), st.normalGradBcCells(), st.normalGradBcCoeffs());

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

    for (const Stencil &st: stencils_)
    {
        //- Boundary assembly
        switch (bType)
        {
            case FIXED:
                eqn.add(st.cell(), st.fixedBcCells(), st.fixedBcCoeffs());
                eqn.addSource(st.cell(), -bRefValue);
                break;

            case ImmersedBoundaryObject::NORMAL_GRADIENT:
                eqn.add(st.cell(), st.normalGradBcCells(), st.normalGradBcCoeffs());
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
            for (const Stencil &st: stencils_)
            {
                eqn.add(st.cell(), st.fixedBcCells(), st.fixedBcCoeffs());
                eqn.addSource(st.cell(), -velocity(st.bp()));
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

//- Force computation

void GhostCellImmersedBoundaryObject::computeForce(Scalar rho,
                                                   Scalar mu,
                                                   const VectorFiniteVolumeField &u,
                                                   const ScalarFiniteVolumeField &p,
                                                   const Vector2D &g)
{
    std::vector<std::tuple<Point2D, Scalar, Scalar>> stresses;
    stresses.reserve(stencils_.size());

    for (const Stencil &st: stencils_)
    {
        stresses.push_back(
                std::make_tuple(
                        st.bp(),
                        st.bpValue(p) + rho * dot(st.bp(), g),
                        mu * dot(dot(st.bpGrad(u), st.nw()), st.nw().tangentVec())
                )
        );
    }

    stresses = grid_->comm().gatherv(grid_->comm().mainProcNo(), stresses);

    if (grid_->comm().isMainProc())
    {
        std::sort(stresses.begin(), stresses.end(),
                  [this](const std::tuple<Point2D, Scalar, Scalar> &a, std::tuple<Point2D, Scalar, Scalar> &b)
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
    {
        stencils_.push_back(Stencil(cell, *this));
        stencils_.back().initFixedBc();
        stencils_.back().initZeroGradientBc();
    }
}
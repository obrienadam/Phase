#include <tuple>

#include "GhostCellImmersedBoundaryObject.h"

GhostCellImmersedBoundaryObject::Stencil::Stencil(const Cell &cell,
                                                  const ImmersedBoundaryObject &ibObj)
        :
        cell_(std::cref(cell))
{
    bp_ = ibObj.nearestIntersect(cell.centroid());
    ip_ = ibObj.isInIb(cell) ? 2. * bp_ - cell.centroid() : 2. * cell.centroid() - bp_;
    nw_ = -ibObj.nearestEdgeNormal(bp_);
}

GhostCellImmersedBoundaryObject::Stencil::Stencil(const Cell &cell,
                                                  const ImmersedBoundaryObject &ibObj,
                                                  const Vector2D &r)
        :
        cell_(std::cref(cell))
{
    bp_ = ibObj.shape().intersections(Ray2D(cell.centroid(), r))[0];
    ip_ = ibObj.isInIb(cell) ? 2. * bp_ - cell.centroid() : 2. * cell.centroid() - bp_;
    nw_ = -ibObj.nearestEdgeNormal(bp_);
}

void GhostCellImmersedBoundaryObject::FixedStencil::init()
{
    bi_.setPoint(ip_);

    if (!bi_.isValid())
        throw Exception("GhostCellImmersedBoundaryObject::FixedStencil", "init", "invalid interpolation stencil.");

    cells_ = bi_.cells();
    cells_.push_back(std::cref(cell_));
    coeffs_ = bi_.coeffs() / 2.;
    coeffs_.push_back(0.5);
}

void GhostCellImmersedBoundaryObject::ZeroGradientStencil::init()
{
    bi_.setPoint(ip_);

    for (const Cell &cell: bi_.cells())
        if (cell_.get().id() == cell.id())
        {
            bi_.setPoint(bp_);
            cells_ = bi_.cells();
            coeffs_ = bi_.derivativeCoeffs(d_);

            if (!bi_.isValid())
                throw Exception("GhostCellImmersedBoundaryObject::ZeroGradientStencil",
                                "init",
                                "invalid interpolation stencil.");

            return;
        }

    if (!bi_.isValid())
        throw Exception("GhostCellImmersedBoundaryObject::ZeroGradientStencil",
                        "init",
                        "invalid interpolation stencil.");

    Scalar l = length();

    cells_ = bi_.cells();
    cells_.push_back(cell_);

    coeffs_ = -bi_.coeffs() / l;
    coeffs_.push_back(1. / l);
}

GhostCellImmersedBoundaryObject::GhostCellImmersedBoundaryObject(const std::string &name,
                                                                 const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                                                 const std::shared_ptr<CellGroup> &solverCells)
        :
        ImmersedBoundaryObject(name, grid, solverCells)
{

}

void GhostCellImmersedBoundaryObject::updateCells()
{
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

    solverCells_->add(cells_);
    cells_.clear();
    ibCells_.clear();
    solidCells_.clear();

    auto items = solverCells_->itemsWithin(*shape_);

    cells_.add(items.begin(), items.end());
    solverCells_->remove(items.begin(), items.end());

    for (const Cell &cell: cells_)
        if (isIbCell(cell))
            ibCells_.add(cell);
        else
            solidCells_.add(cell);

    constructStencils();
}

FiniteVolumeEquation<Scalar> GhostCellImmersedBoundaryObject::bcs(ScalarFiniteVolumeField &field) const
{
    FiniteVolumeEquation<Scalar> eqn(field);
    BoundaryType bType = boundaryType(field.name());
    Scalar bRefValue = getBoundaryRefValue<Scalar>(field.name());

    switch (bType)
    {
        case FIXED:
            for (const auto &st: fixedStencils_)
            {
                eqn.add(st.cell(), st.cells(), st.coeffs());
                eqn.addSource(st.cell(), -bRefValue);
            }

            for (const Cell &cell: solidCells_)
            {
                eqn.add(cell, cell, 1.);
                eqn.addSource(cell, -bRefValue);
            }

            break;
        case NORMAL_GRADIENT:
            for (const Stencil &st: zeroGradientStencils_)
                eqn.add(st.cell(), st.cells(), st.coeffs());

            for (const Cell &cell: solidCells_)
                eqn.add(cell, cell, 1.);
            break;

        default:
            throw Exception("GhostCellImmersedBoundaryObject", "bcs", "invalid boundary type.");
    }

    return eqn;
}

FiniteVolumeEquation<Vector2D> GhostCellImmersedBoundaryObject::bcs(VectorFiniteVolumeField &field) const
{
    FiniteVolumeEquation<Vector2D> eqn(field);
    BoundaryType bType = boundaryType(field.name());
    Scalar bRefValue = getBoundaryRefValue<Scalar>(field.name());

    //- Boundary assembly
    switch (bType)
    {
        case FIXED:
            for (const auto &st: fixedStencils_)
            {
                eqn.add(st.cell(), st.cells(), st.coeffs());
                eqn.addSource(st.cell(), -bRefValue);
            }
            break;

        case ImmersedBoundaryObject::NORMAL_GRADIENT:
            for (const Stencil &st: zeroGradientStencils_)
            {
                eqn.add(st.cell(), st.cells(), st.coeffs());
                eqn.addSource(st.cell(), -bRefValue);
            }
            break;

        default:
            throw Exception("GhostCellImmersedBoundaryObject", "bcs", "invalid boundary type.");
    }

    for (const Cell &cell: solidCells_)
    {
        eqn.add(cell, cell, 1.);
        eqn.addSource(cell, -bRefValue);
    }

    return eqn;
}

FiniteVolumeEquation<Vector2D> GhostCellImmersedBoundaryObject::velocityBcs(VectorFiniteVolumeField &u) const
{
    FiniteVolumeEquation<Vector2D> eqn(u);

    switch (boundaryType(u.name()))
    {
        case FIXED:
            for (const auto &st: fixedStencils_)
            {
                eqn.add(st.cell(), st.cells(), st.coeffs());
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

FiniteVolumeEquation<Scalar> GhostCellImmersedBoundaryObject::pressureBcs(ScalarFiniteVolumeField &p) const
{
    FiniteVolumeEquation<Scalar> eqn(p);

    for (const auto &st: zeroGradientStencils_)
    {
        eqn.add(st.cell(), st.cells(), st.coeffs());
        eqn.addSource(st.cell(), dot(acceleration(st.bp()), st.nw()));
    }

    for (const Cell &cell: solidCells_)
    {
        eqn.add(cell, cell, 1.);
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
    stresses.reserve(fixedStencils_.size());

    for (const auto &st: fixedStencils_)
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
    fixedStencils_.clear();
    zeroGradientStencils_.clear();

    for (const Cell &cell: ibCells_)
    {
        fixedStencils_.push_back(FixedStencil(cell, *this));
        fixedStencils_.back().init();

        zeroGradientStencils_.push_back(ZeroGradientStencil(cell, *this));
        zeroGradientStencils_.back().init();
    }
}
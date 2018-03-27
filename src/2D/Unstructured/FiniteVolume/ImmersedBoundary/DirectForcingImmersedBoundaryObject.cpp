#include "Math/Matrix.h"
#include "FiniteVolumeGrid2D/BilinearInterpolator.h"

#include "DirectForcingImmersedBoundaryObject.h"

DirectForcingImmersedBoundaryObject::Stencil::Stencil(const VectorFiniteVolumeField &u,
                                                      const Cell &cell,
                                                      const ImmersedBoundaryObject &ibObj)
{
    bp_ = ibObj.nearestIntersect(cell.centroid());
    ip_ = 2. * cell.centroid() - bp_;

    ub_ = ibObj.velocity(bp_);
    uip_ = BilinearInterpolator(ibObj.grid(), ip_)(u);
    uf_ = (ub_ + uip_) / 2.;
}

DirectForcingImmersedBoundaryObject::FieldExtensionStencil::FieldExtensionStencil(const VectorFiniteVolumeField &u,
                                                                                  const Cell &cell,
                                                                                  const ImmersedBoundaryObject &ibObj)
{
    bp_ = ibObj.nearestIntersect(cell.centroid());
    ip_ = 2. * bp_ - cell.centroid();

    BilinearInterpolator bi(ibObj.grid(), ip_);

    uip_ = bi(u);
    ub_ = ibObj.velocity(bp_);
    uf_ = 2. * ub_ - uip_;

    cells_ = bi.cells();
    cells_.push_back(std::cref(cell));

    coeffs_ = bi.coeffs() / 2.;
    coeffs_.push_back(0.5);
}


DirectForcingImmersedBoundaryObject::PressureFieldExtensionStencil::PressureFieldExtensionStencil(
        const ScalarFiniteVolumeField &p,
        const Cell &cell,
        const ImmersedBoundaryObject &ibObj)
{
    bp_ = ibObj.nearestIntersect(cell.centroid());
    ip_ = 2. * bp_ - cell.centroid();
    n_ = ibObj.nearestEdgeNormal(bp_).unitVec();

    BilinearInterpolator bi(ibObj.grid(), bp_);

    bool validStencil = false;

    for(const Cell& kCell: cells_)
        if(kCell.id() == cell.id())
        {
            validStencil = true;
            break;
        }

    if(validStencil)
    {
        cells_ = bi.cells();
        coeffs_ = bi.derivativeCoeffs(n_);
    }
    else
    {
        bi.setPoint(ip_);
        n_ = (cell.centroid() - ip_).unitVec();
        cells_ = bi.cells();
        cells_.push_back(std::cref(cell));
        coeffs_ = -bi.coeffs() / (cell.centroid() - ip_).mag();
        coeffs_.push_back(1. / (cell.centroid() - ip_).mag());
    }
}

DirectForcingImmersedBoundaryObject::QuadraticStencil::QuadraticStencil(const VectorFiniteVolumeField &u,
                                                                        const Cell &cell,
                                                                        const ImmersedBoundaryObject &ibObj)
{
    std::vector<Ref<const Cell>> cells;
    std::vector<Point2D> bps;

    auto isIbCell = [&ibObj](const Cell &cell)
    {
        for (const InteriorLink &nb: cell.neighbours())
            if (ibObj.isInIb(nb.cell()))
                return true;

        return false;
    };

    for (const CellLink &nb: cell.cellLinks())
    {
        if (ibObj.isInIb(nb.cell()))
            continue;

        if (isIbCell(nb.cell()))
        {
            cells.push_back(std::cref(nb.cell()));
            bps.push_back(ibObj.nearestIntersect(nb.cell().centroid()));
        }
        else
            cells.push_back(std::cref(nb.cell()));
    }

    bps.push_back(ibObj.nearestIntersect(cell.centroid()));

    Matrix A(cells.size() + bps.size(), 6);
    Matrix b(cells.size() + bps.size(), 2);

    int i = 0;
    for (const Cell &cell: cells)
    {
        Point2D x = cell.centroid();
        A(i, 0) = x.x * x.x;
        A(i, 1) = x.y * x.y;
        A(i, 2) = x.x * x.y;
        A(i, 3) = x.x;
        A(i, 4) = x.y;
        A(i, 5) = 1.;
        b(i, 0) = u(cell).x;
        b(i++, 1) = u(cell).y;
    }

    for (const Point2D &x: bps)
    {
        A(i, 0) = x.x * x.x;
        A(i, 1) = x.y * x.y;
        A(i, 2) = x.x * x.y;
        A(i, 3) = x.x;
        A(i, 4) = x.y;
        A(i, 5) = 1.;

        Vector2D ub = ibObj.velocity(x);
        b(i, 0) = ub.x;
        b(i++, 1) = ub.y;
    }

    Point2D x = cell.centroid();
    Matrix tmp = Matrix(1, 6, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.}) * solve(A, b);
    uf_ = Vector2D(tmp(0, 0), tmp(0, 1));
}

DirectForcingImmersedBoundaryObject::DirectForcingImmersedBoundaryObject(const std::string &name,
                                                                         Label id,
                                                                         const std::shared_ptr<FiniteVolumeGrid2D> &grid)
        :
        ImmersedBoundaryObject(name, id, grid)
{

}

void DirectForcingImmersedBoundaryObject::updateCells()
{
    ibCells_.clear();
    solidCells_.clear();
    forcingCells_.clear();
    pseudoFluidCells_.clear();

    // auto items = fluid_->itemsWithin(*shape_);

    for (const Cell &cell: fluid_->itemsWithin(*shape_))
    {
        solidCells_.add(cell);

        for (const CellLink &nb: cell.neighbours())
        {
            if (!shape_->isInside(nb.cell().centroid()))
            {
                ibCells_.add(nb.cell());
                forcingCells_.add(nb.cell());
                pseudoFluidCells_.add(cell);
            }
        }
    }
}

void DirectForcingImmersedBoundaryObject::computeBoundaryForcing(const VectorFiniteVolumeField &u,
                                                                 Scalar timeStep,
                                                                 VectorFiniteVolumeField &fb) const
{
    for (const Cell &cell: ibCells_)
    {
        Stencil st = Stencil(u, cell, *this);
        fb(cell) = (st.uf() - u(cell)) / timeStep;
    }

    for (const Cell &cell: solidCells_)
    {
        fb(cell) = (velocity(cell.centroid()) - u(cell)) / timeStep;
    }
}
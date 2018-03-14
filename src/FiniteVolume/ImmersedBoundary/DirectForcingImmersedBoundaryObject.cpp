#include "DirectForcingImmersedBoundaryObject.h"

DirectForcingImmersedBoundaryObject::Stencil::Stencil(const VectorFiniteVolumeField &u,
                                                      const Cell &cell,
                                                      const ImmersedBoundaryObject &ibObj)
{
    bp_ = ibObj.nearestIntersect(cell.centroid());
    ip_ = 2. * cell.centroid() - bp_;
    ipCells_ = ibObj.grid()->findNearestNode(ip_).cells();
    ub_ = ibObj.velocity(bp_);

    Point2D x[] = {
            ipCells_[0].get().centroid(),
            ipCells_[1].get().centroid(),
            ipCells_[2].get().centroid(),
            ipCells_[3].get().centroid()
    };

    StaticMatrix<4, 4> A = {
            x[0].x * x[0].y, x[0].x, x[0].y, 1.,
            x[1].x * x[1].y, x[1].x, x[1].y, 1.,
            x[2].x * x[2].y, x[2].x, x[2].y, 1.,
            x[3].x * x[3].y, x[3].x, x[3].y, 1.
    };

    StaticMatrix<4, 2> b = {
            u(ipCells_[0]).x, u(ipCells_[0]).y,
            u(ipCells_[1]).x, u(ipCells_[1]).y,
            u(ipCells_[2]).x, u(ipCells_[2]).y,
            u(ipCells_[3]).x, u(ipCells_[3]).y,
    };

    StaticMatrix<1, 2> c = StaticMatrix<1, 4>({ip_.x * ip_.y, ip_.x, ip_.y, 1.}) * solve(A, b);

    uip_ = Vector2D(c(0, 0), c(0, 1));
    uf_ = (ub_ + uip_) / 2.;

    StaticMatrix<1, 4> a = StaticMatrix<1, 4>({ip_.x * ip_.y, ip_.x, ip_.y, 1.}) * inverse(A);
    ipCoeffs_.assign(a.begin(), a.end());
}

DirectForcingImmersedBoundaryObject::FieldExtensionStencil::FieldExtensionStencil(const VectorFiniteVolumeField &u,
                                                                                  const Cell &cell,
                                                                                  const ImmersedBoundaryObject &ibObj)
{
    bp_ = ibObj.nearestIntersect(cell.centroid());
    ip_ = 2. * bp_ - cell.centroid();
    ipCells_ = ibObj.grid()->findNearestNode(ip_).cells();

    Point2D x[] = {
            ipCells_[0].get().centroid(),
            ipCells_[1].get().centroid(),
            ipCells_[2].get().centroid(),
            ipCells_[3].get().centroid()
    };

    StaticMatrix<4, 4> A = {
            x[0].x * x[0].y, x[0].x, x[0].y, 1.,
            x[1].x * x[1].y, x[1].x, x[1].y, 1.,
            x[2].x * x[2].y, x[2].x, x[2].y, 1.,
            x[3].x * x[3].y, x[3].x, x[3].y, 1.
    };

    StaticMatrix<4, 2> b = {
            u(ipCells_[0]).x, u(ipCells_[0]).y,
            u(ipCells_[1]).x, u(ipCells_[1]).y,
            u(ipCells_[2]).x, u(ipCells_[2]).y,
            u(ipCells_[3]).x, u(ipCells_[3]).y,
    };

    StaticMatrix<1, 2> c = StaticMatrix<1, 4>({ip_.x * ip_.y, ip_.x, ip_.y, 1.}) * solve(A, b);

    uip_ = Vector2D(c(0, 0), c(0, 1));
    ub_ = ibObj.velocity(bp_);
    uf_ = 2. * ub_ - uip_;

    StaticMatrix<1, 4> a = StaticMatrix<1, 4>({ip_.x * ip_.y, ip_.x, ip_.y, 1.}) * inverse(A);
    ipCoeffs_.assign(a.begin(), a.end());
}

DirectForcingImmersedBoundaryObject::QuadraticStencil::QuadraticStencil(const VectorFiniteVolumeField &u,
                                                                        const Cell &cell,
                                                                        const ImmersedBoundaryObject &ibObj)
{
    std::vector<Ref<const Cell>> cells;
    std::vector<Point2D> bps;

    auto isIbCell = [&ibObj](const Cell &cell) {
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

DirectForcingImmersedBoundaryObject::DirectForcingImmersedBoundaryObject(const std::string &name, Label id,
                                                                         const ImmersedBoundary &ib,
                                                                         const std::shared_ptr<FiniteVolumeGrid2D> &grid)
        :
        ImmersedBoundaryObject(name, id, ib, grid)
{

}

void DirectForcingImmersedBoundaryObject::updateCells()
{
    ibCells_.clear();
    solidCells_.clear();

    for (const Cell &cell: fluid_->itemsWithin(*shape_))
    {
        for (const CellLink &nb: cell.neighbours())
        {
            if (!shape_->isInside(nb.cell().centroid()))
            {
                ibCells_.add(nb.cell());
                solidCells_.add(cell);
            }
        }
    }
}

void DirectForcingImmersedBoundaryObject::updateIbForce(const VectorFiniteVolumeField &u,
                                                        Scalar timeStep,
                                                        VectorFiniteVolumeField &fb)
{
    fb.fill(Vector2D(0., 0.));

    for (const Cell &cell: ibCells_)
    {
        Stencil st = Stencil(u, cell, *this);
        fb(cell) = (st.uf() - u(cell)) / timeStep;
    }
}
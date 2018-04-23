#include "Math/Matrix.h"
#include "FiniteVolumeGrid2D/BilinearInterpolator.h"

#include "DirectForcingImmersedBoundaryObject.h"

DirectForcingImmersedBoundaryObject::Stencil::Stencil(const VectorFiniteVolumeField &u,
                                                      const Cell &cell,
                                                      const ImmersedBoundaryObject &ibObj)
{
    bp_ = ibObj.nearestIntersect(cell.centroid());
    ub_ = ibObj.velocity(bp_);

    auto isIbCell = [&ibObj](const Cell& cell)
    {
        for(const CellLink &nb: cell.neighbours())
            if(ibObj.isInIb(nb.cell()))
                return true;
        return false;
    };

    std::vector<Ref<const Cell>> iCells;

    for(const CellLink &nb: cell.cellLinks())
        if(!ibObj.isInIb(nb.cell()) && !isIbCell(nb.cell()))
        {
            iCells.push_back(std::cref(nb.cell()));
        }

    if(iCells.size() < 2)
    {
        throw Exception("DirectForcingImmersedBoundaryObject::Stencil", "Stencil",
                        "not enough cells to form stencil. N cells = " + std::to_string(cell.cellLinks().size())
                        + ", proc = " + std::to_string(u.grid()->comm().rank())
                        + ", cell id = " + std::to_string(u.grid()->globalIds()[cell.id()]) + ".");
    }

    Point2D x[] = {iCells[0].get().centroid(), iCells[1].get().centroid()};
    Vector2D iu[] = {u(iCells[0]), u(iCells[1])};

    auto c = solve<3, 2>({
                             x[0].x, x[0].y, 1.,
                             x[1].x, x[1].y, 1.,
                             bp_.x, bp_.y, 1.
                         }, {
                             iu[0].x, iu[0].y,
                             iu[1].x, iu[1].y,
                             ub_.x, ub_.y
                         });

    Point2D xc = cell.centroid();

    uf_ = Vector2D(c(0, 0) * xc.x + c(1, 0) * xc.y + c(2, 0),
                   c(0, 1) * xc.x + c(1, 1) * xc.y + c(2, 1));
}

DirectForcingImmersedBoundaryObject::QuadraticLsStencil::QuadraticLsStencil(const VectorFiniteVolumeField &u, const Cell &cell, const ImmersedBoundaryObject &ibObj)
{
    bp_ = ibObj.nearestIntersect(cell.centroid());
    ub_ = ibObj.velocity(cell.centroid());

    std::vector<Point2D> x;
    std::vector<Vector2D> src;

    for(const CellLink &nb: cell.cellLinks())
        if(!ibObj.isInIb(nb.cell()))
        {
            x.push_back(nb.cell().centroid());
            src.push_back(u(nb.cell()));

            for(const CellLink &nb2: nb.cell().neighbours())
            {
                if(ibObj.isInIb(nb2.cell().centroid()))
                {
                    Point2D compat = ibObj.nearestIntersect(nb.cell().centroid());
                    x.push_back(compat);
                    src.push_back(ibObj.velocity(compat));
                    break;
                }
            }
        }

    x.push_back(bp_);
    src.push_back(ub_);

    Matrix A(x.size(), 6), b(x.size(), 2);

    for(int i = 0; i < x.size(); ++i)
    {
        A.setRow(i, {x[i].x * x[i].x, x[i].y * x[i].y, x[i].x * x[i].y, x[i].x, x[i].y, 1.});
        b.setRow(i, {src[i].x, src[i].y});
    }

    Point2D xc = cell.centroid();

    Matrix sln = Matrix(1, 6, {xc.x * xc.x, xc.y * xc.y, xc.x * xc.y, xc.x, xc.y, 1.}) * solve(A, b);
    uf_ = Vector2D(sln(0, 0), sln(0, 1));
}

DirectForcingImmersedBoundaryObject::DirectForcingImmersedBoundaryObject(const std::string &name,
                                                                         const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                                                         const std::shared_ptr<CellGroup> &solverCells)
    :
      ImmersedBoundaryObject(name, grid, solverCells)
{

}

void DirectForcingImmersedBoundaryObject::updateCells()
{
    solidCells_.clear();
    ibCells_.clear();

    auto items = solverCells_->itemsWithin(*shape_);

    solidCells_.add(items.begin(), items.end());

    for(const Cell& cell: solidCells_)
        for(const CellLink& nb: cell.neighbours())
            if(!isInIb(nb.cell()) && grid_->localCells().isInGroup(nb.cell()))
                ibCells_.add(nb.cell());
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

FiniteVolumeEquation<Vector2D> DirectForcingImmersedBoundaryObject::velocityBcs(VectorFiniteVolumeField &u) const
{
    FiniteVolumeEquation<Vector2D> eqn(u);

    Scalar timeStep = 5e-3;

    for(const Cell& cell: ibCells_)
    {

        std::vector<Ref<const Cell>> iCells;
        std::vector<Point2D> bPts;
        std::vector<Vector2D> bVel;

        for(const CellLink &nb: cell.cellLinks())
            if(!isInIb(nb.cell()))
            {
                iCells.push_back(std::cref(nb.cell()));

                for(const CellLink &nb2: nb.cell().neighbours())
                {
                    if(isInIb(nb2.cell().centroid()))
                    {
                        bPts.push_back(nearestIntersect(nb.cell().centroid()));
                        bVel.push_back(velocity(bPts.back()));
                        break;
                    }
                }
            }

        bPts.push_back(nearestIntersect(cell.centroid()));
        bVel.push_back(velocity(bPts.back()));

        Matrix A(iCells.size() + bPts.size(), 6);

        for(int i = 0; i < iCells.size(); ++i)
        {
            Point2D x = iCells[i].get().centroid();
            A.setRow(i, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});
        }

        for(int i = 0; i < bPts.size(); ++i)
        {
            Point2D x = bPts[i];
            A.setRow(i + iCells.size(), {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});
        }

        Point2D xc = cell.centroid();
        Matrix x = Matrix(1, 6, {xc.x * xc.x, xc.y * xc.y, xc.x * xc.y, xc.x, xc.y, 1.}) * pseudoInverse(A) * cell.volume() / timeStep;

        eqn.add(cell, cell, -cell.volume() / timeStep);

        for(int i = 0; i < iCells.size(); ++i)
        {
            eqn.add(cell, iCells[i], x(0, i));
        }

        for(int i = 0; i < bPts.size(); ++i)
        {
            eqn.addSource(cell, x(0, i + iCells.size()) * bVel[i]);
        }
    }

    for(const Cell& cell: solidCells_)
    {
        eqn.add(cell, cell, -cell.volume() / timeStep);
        eqn.addSource(cell, cell.volume() * velocity(cell.centroid()) / timeStep);
    }

    return eqn;
}

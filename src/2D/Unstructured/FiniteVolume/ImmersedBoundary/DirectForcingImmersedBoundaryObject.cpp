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

DirectForcingImmersedBoundaryObject::FieldExtensionStencil::FieldExtensionStencil(const Cell &cell, const ImmersedBoundaryObject &ibObj)
{
    cell_ = &cell;
    bp_ = ibObj.nearestIntersect(cell.centroid());
    nb_ = ibObj.nearestEdgeNormal(bp_).normalVec();
    ub_ = ibObj.velocity(bp_);
    ab_ = ibObj.acceleration(bp_);
    iCells_.reserve(5);

    for(const CellLink &nb: cell.cellLinks())
        if(!ibObj.isInIb(nb.cell()))
            iCells_.push_back(&nb.cell());

    if(iCells_.size() < 2)
    {
        throw Exception("DirectForcingImmersedBoundaryObject::FieldExtensionStencil", "FieldExtensionStencil",
                        "not enough cells to form stencil. N cells = " + std::to_string(cell.cellLinks().size())
                        + ", proc = " + std::to_string(ibObj.grid()->comm().rank())
                        + ", cell id = " + std::to_string(ibObj.grid()->globalIds()[cell.id()]) + ".");
    }

    Point2D x[] = {iCells_[0]->centroid(), iCells_[1]->centroid()};

    Au_ = inverse<3>({
                         x[0].x, x[0].y, 1.,
                         x[1].x, x[1].y, 1.,
                         bp_.x, bp_.y, 1.
                     });

    Ap_ = inverse<3>({
                         x[0].x, x[0].y, 1.,
                         x[1].x, x[1].y, 1.,
                         nb_.x, nb_.y, 0.
                     });
}

Vector2D DirectForcingImmersedBoundaryObject::FieldExtensionStencil::uExtend(const VectorFiniteVolumeField &u) const
{
    auto c = Au_ * StaticMatrix<3, 2>({
                                          u(*iCells_[0]).x, u(*iCells_[0]).y,
                                          u(*iCells_[1]).x, u(*iCells_[1]).y,
                                          ub_.x, ub_.y
                                      });

    return Vector2D(
                c(0, 0) * cell_->centroid().x + c(1, 0) * cell_->centroid().y + c(2, 0),
                c(0, 1) * cell_->centroid().x + c(1, 1) * cell_->centroid().y + c(2, 1)
                );
}

Scalar DirectForcingImmersedBoundaryObject::FieldExtensionStencil::pExtend(const ScalarFiniteVolumeField &p) const
{
    auto c = Ap_ * StaticMatrix<3, 1>({
                                          p(*iCells_[0]),
                                          p(*iCells_[1]),
                                          -dot(ab_, nb_)
                                      });

    return c(0, 0) * cell_->centroid().x + c(1, 0) * cell_->centroid().y + c(2, 0);
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
    return eqn;
}

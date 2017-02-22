#include "OscillatingImmersedBoundaryObject.h"

OscillatingImmersedBoundaryObject::OscillatingImmersedBoundaryObject(const std::string &name,
                                                                     Scalar amplitude,
                                                                     Scalar frequency,
                                                                     Label id,
                                                                     FiniteVolumeGrid2D &grid)
    :
      ImmersedBoundaryObject(name, id, grid)
{
    amp_ = amplitude;
    omega_ = 2.*M_PI*frequency;
}

void OscillatingImmersedBoundaryObject::update(Scalar timeStep)
{
    if(time_ == 0.)
        origin_ = shape().centroid();

    Vector2D pos = origin_ + amp_*sin(omega_*(time_ += timeStep));
    shape() += pos - shape().centroid();
    updateCells();
}

void OscillatingImmersedBoundaryObject::updateCells()
{
    CellZone& fluidCells = grid_.cellZone("fluid");

    for(Cell& cell: *freshlyClearedCells_)
        fluidCells.moveToGroup(cell);

    grid_.setCellsActive(*freshlyClearedCells_);
    freshlyClearedCells_->clear();

    for(Cell& cell: *cells_)
        if(!shapePtr_->isBoundedBy(cell.centroid(), 1e-10))
            freshlyClearedCells_->push_back(cell);

    ibCells_->clear();
    solidCells_->clear();

    for(Cell& cell: fluidCells.rangeSearch(*shapePtr_, 1e-10))
        cells_->moveToGroup(cell);

    auto isIbCell = [this](const Cell& cell) {
        for(const InteriorLink& nb: cell.neighbours())
            if(!shapePtr_->isBoundedBy(nb.cell().centroid(), 1e-10))
                return true;

        for(const DiagonalCellLink& dg: cell.diagonals())
            if(!shapePtr_->isBoundedBy(dg.cell().centroid(), 1e-10))
                return true;

        return false;
    };

    for(Cell& cell: *cells_)
        if(isIbCell(cell))
            ibCells_->push_back(cell);
        else
            solidCells_->push_back(cell);

    grid_.setCellsActive(*ibCells_);
    grid_.setCellsInactive(*solidCells_);
    constructStencils();
}

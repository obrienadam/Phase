#include "TranslatingImmersedBoundaryObject.h"


TranslatingImmersedBoundaryObject::TranslatingImmersedBoundaryObject(const std::string &name,
                                                                     const Vector2D &velocity,
                                                                     Scalar omega,
                                                                     Label id,
                                                                     FiniteVolumeGrid2D &grid)
    :
      ImmersedBoundaryObject(name,
                             id,
                             grid),
      velocity_(velocity),
      omega_(omega)
{

}

void TranslatingImmersedBoundaryObject::update(Scalar timeStep)
{
    shape() += timeStep*velocity_;
    updateCells();
}

void TranslatingImmersedBoundaryObject::updateCells()
{
    CellZone& fluidCells = grid_.cellZone("fluid");
    freshlyClearedCells_->clear();

    for(Cell& cell: *cells_)
        if(!shapePtr_->isBoundedBy(cell.centroid(), 1e-10))
        {
            freshlyClearedCells_->push_back(cell); //- Freshly cleared cells need a correction
            fluidCells.moveToGroup(cell); //- This will remove the cell from the ib cell zone
        }

    grid_.setCellsActive(*freshlyClearedCells_);

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

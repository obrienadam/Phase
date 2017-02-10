#include "MovingImmersedBoundaryObject.h"

void MovingImmersedBoundaryObject::updateCells()
{
    CellZone& fluidZone = grid_.cellZone("fluid");
    std::vector<Ref<const Cell>> cells = cells_->cells();

    grid_.setCellsActive(grid_.getCellIds(cells));

    for(const Cell& cell: cells)
        fluidZone.moveToGroup(cell);

    cells_->clear();
    ibCells_->clear();
    solidCells_->clear();

    flagIbCells();
    constructStencils();
}

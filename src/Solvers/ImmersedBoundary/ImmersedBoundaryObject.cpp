#include "ImmersedBoundaryObject.h"
#include "CellSearch.h"

//- Immersed boundary object stencil class

ImmersedBoundaryObject::ImmersedBoundaryStencil::ImmersedBoundaryStencil(const Cell &cell,
                                                                         const Point2D &bp,
                                                                         const FiniteVolumeGrid2D &grid,
                                                                         const CellSearch &cs)
    :
      cell_(cell),
      bp_(bp),
      ip_(2.*bp - cell.centroid())
{
    Point2D sp = grid.findNearestNode(ip_);
    kNN_ = cs.kNearestNeighbourSearch(sp, 4);
}

//- Immersed boundary objectclass

ImmersedBoundaryObject::ImmersedBoundaryObject(const FiniteVolumeGrid2D &grid)
    :
      Circle(Point2D(), 0.),
      grid_(grid)
{

}

ImmersedBoundaryObject::ImmersedBoundaryObject(const FiniteVolumeGrid2D &grid, const Point2D &center, Scalar radius)
    :
      Circle(center, radius),
      grid_(grid)
{

}

void ImmersedBoundaryObject::init(const Point2D& center, Scalar radius)
{
    Circle::init(center, radius);
}

void ImmersedBoundaryObject::constructStencils()
{
    grid().moveAllCellsToFluidCellGroup(); // This must happen first for some reason

    // Get a list of cells that are in immersed boundary
    std::vector< Ref<const Cell> > internalCells = CellSearch(grid_.activeCells()).rangeSearch(*this);

    ibStencils_.clear();
    ibStencils_.reserve(internalCells.size());

    std::vector<size_t> inactiveCellIds, ibCellIds;

    inactiveCellIds.reserve(internalCells.size());
    ibCellIds.reserve(internalCells.size());

    for(const Cell &cell: internalCells)
    {
        bool isIbCell = false;

        for(const InteriorLink &nb: cell.neighbours())
        {
            if(!isInside(nb.cell().centroid()))
            {
                isIbCell = true;
                break;
            }
        }

        if(isIbCell)
        {
            ibCellIds.push_back(cell.id()); // set the cell as an ib cell (part of a new group)
            ibStencils_.push_back(ImmersedBoundaryStencil(
                                      cell,
                                      nearestIntersect(cell.centroid()),
                                      grid_,
                                      grid_.activeCells()));
        }
        else // the cell will be set to inactive
            inactiveCellIds.push_back(cell.id());
    }

    grid_.moveCellsToInactiveCellGroup(inactiveCellIds);
    grid_.moveCellsToCellGroup("ibCells", ibCellIds);
}

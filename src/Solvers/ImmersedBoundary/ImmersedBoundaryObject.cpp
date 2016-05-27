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
      ip_(2.*bp - cell.centroid()),
      cellSearch_(cs)
{
    Point2D sp = grid.findNearestNode(ip_);
    kNN_ = cs.kNearestNeighbourSearch(sp, 4);
}

ImmersedBoundaryObject::ImmersedBoundaryStencil ImmersedBoundaryObject::ImmersedBoundaryStencil::rotate(Scalar theta, const ImmersedBoundaryObject& ibObj) const
{
    Vector2D bp = bp_ - ibObj.centroid();
    bp = bp.rotate(theta) + ibObj.centroid();

    return ImmersedBoundaryStencil(cell_, bp, ibObj.grid(), cellSearch_);
}

//- Immersed boundary objectclass

ImmersedBoundaryObject::ImmersedBoundaryObject(const FiniteVolumeGrid2D &grid, const std::shared_ptr<SurfaceTensionForce> &csfPtr)
    :
      Circle(Point2D(), 0.),
      grid_(grid),
      cellSearch_(grid.activeCells())
{
    csf_ = csfPtr;
}

ImmersedBoundaryObject::ImmersedBoundaryObject(const FiniteVolumeGrid2D &grid, const std::shared_ptr<SurfaceTensionForce> &csfPtr, const Point2D &center, Scalar radius)
    :
      Circle(center, radius),
      grid_(grid),
      cellSearch_(grid.activeCells())
{
    csf_ = csfPtr;
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
                                      cellSearch_
                                      )
                                  );
        }
        else // the cell will be set to inactive
            inactiveCellIds.push_back(cell.id());
    }

    grid_.moveCellsToInactiveCellGroup(inactiveCellIds);
    grid_.moveCellsToCellGroup("ibCells", ibCellIds);
}

void ImmersedBoundaryObject::addBoundaryType(const std::string &fieldName, BoundaryType boundaryType)
{
    boundaryTypes_[fieldName] = boundaryType;
}

std::vector< Ref<const Cell> > ImmersedBoundaryObject::getBoundingCells(const Point2D &pt) const
{
    Vector2D node = grid_.findNearestNode(pt);
    return cellSearch_.kNearestNeighbourSearch(node, 4);
}

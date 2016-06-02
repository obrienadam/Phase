#include "ImmersedBoundaryObject.h"

ImmersedBoundaryObject::ImmersedBoundaryObject(const FiniteVolumeGrid2D& grid, const ImmersedBoundaryContinuumSurfaceForce &csf)
    :
      Circle(Point2D(0., 0.), 0.),
      grid_(grid),
      csf_(csf)
{

}

Point2D ImmersedBoundaryObject::boundaryPoint(const Point2D &pt) const
{
    return nearestIntersect(pt);
}

Point2D ImmersedBoundaryObject::imagePoint(const Point2D &pt) const
{
    return 2.*boundaryPoint(pt) - pt;
}

void ImmersedBoundaryObject::setInternalCells()
{
    auto internalCells = grid_.activeCells().rangeSearch(*this);
    std::vector<size_t> internalCellIds;
    internalCellIds.reserve(internalCells.size());

    grid_.moveAllCellsToFluidCellGroup(); // temporary!

    for(const Cell &cell: internalCells)
        internalCellIds.push_back(cell.id());

    grid_.moveCellsToInactiveCellGroup(internalCellIds);
    internalCellIds.clear();

    for(const Cell &cell: internalCells)
    {
        for(const InteriorLink &nb: cell.neighbours())
        {
            if(nb.cell().isActive())
            {
                internalCellIds.push_back(cell.id());
                break;
            }
        }
    }

    grid_.moveCellsToCellGroup("ibCells", internalCellIds);
}

void ImmersedBoundaryObject::addBoundaryType(const std::string &name, BoundaryType boundaryType)
{
    boundaryTypes_[name] = boundaryType;
}

const std::vector< Ref<const Cell> > ImmersedBoundaryObject::boundingCells(const Point2D &pt) const
{
    const Node &node = grid_.findNearestNode(pt);
    return grid_.activeCells().kNearestNeighbourSearch(node, 4);
}

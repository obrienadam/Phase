#include "ImmersedBoundaryObject.h"

ImmersedBoundaryObject::ImmersedBoundaryObject(const std::string& name,
                                               const FiniteVolumeGrid2D& grid,
                                               const Point2D& center,
                                               Scalar radius)
    :
      Circle(center, radius),
      name_(name),
      grid_(grid)
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

    for(const Cell &cell: internalCells)
        internalCellIds.push_back(cell.id());

    grid_.moveCellsToInactiveCellGroup(internalCellIds);
    internalCellIds.clear();

    for(const Cell &cell: internalCells)
    {
        bool isActive = false;

        for(const InteriorLink &nb: cell.neighbours())
        {
            if(nb.cell().isActive())
            {
                internalCellIds.push_back(cell.id());
                isActive = true;
                break;
            }
        }

        if(isActive)
            continue;

        for(const DiagonalCellLink &dg: cell.diagonals())
        {
            if(dg.cell().isActive())
            {
                internalCellIds.push_back(cell.id());
                break;
            }
        }
    }

    grid_.moveCellsToCellGroup(name_ + "_cells", internalCellIds);
}

void ImmersedBoundaryObject::addBoundaryType(const std::string &name, BoundaryType boundaryType)
{
    boundaryTypes_[name] = boundaryType;
}

void ImmersedBoundaryObject::addBoundaryRefValue(const std::string &name, Scalar boundaryRefValue)
{
    boundaryRefValues_[name] = boundaryRefValue;
}

const std::vector< Ref<const Cell> > ImmersedBoundaryObject::boundingCells(const Point2D &pt) const
{
    const Node &node = grid_.findNearestNode(pt);
    return grid_.activeCells().kNearestNeighbourSearch(node, 4);
}

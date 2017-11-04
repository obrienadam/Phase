#include "CellLink.h"
#include "Cell.h"

CellLink::CellLink(const Cell &self, const Cell &other)
        :
        Link(self),
        cell_(other)
{
    rCellVec_ = cell_.centroid() - self_.centroid();
}
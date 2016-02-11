#include <string>

#include "FiniteVolumeGrid2D.h"
#include "Exception.h"

FiniteVolumeGrid2D::FiniteVolumeGrid2D(int nCellsI, int nCellsJ)
    :
      cornerNodes_(nCellsI + 1, nCellsJ + 1, "CornerNodes"),
      cellNodes_(nCellsI, nCellsJ, 0., "CellNodes"),
      faceNodesI_(nCellsI + 1, nCellsJ, "FaceNodesI"),
      faceNodesJ_(nCellsI, nCellsJ + 1, "FaceNodesJ"),
      cellIndices_(nCellsI, nCellsJ, "CellIndices")
{
    initCellIndices();
}

const Point2D& FiniteVolumeGrid2D::cornerNode(int i, int j, Node node) const
{
    switch (node)
    {
    case SW: return cornerNodes_(i, j);
    case SE: return cornerNodes_(i + 1, j);
    case NE: return cornerNodes_(i + 1, j + 1);
    case NW: return cornerNodes_(i, j + 1);
    default: throw Exception("FiniteVolumeGrid2D", "cornerNode", "invalid node.");
    };
}

bool FiniteVolumeGrid2D::inRange(int i, int j) const
{
    return i >= 0 && j >= 0 && i < nCellsI() && j < nCellsJ();
}

int FiniteVolumeGrid2D::cellIndex(int i, int j) const
{
    return cellIndices_(i, j);
}

void FiniteVolumeGrid2D::initCellIndices()
{
    int id = 0;
    nActiveCells_ = 0;

    for(int j = 0; j < cellIndices_.sizeJ(); ++j)
        for(int i = 0; i < cellIndices_.sizeI(); ++i)
        {
            cellIndices_(i, j) = id++;
            ++nActiveCells_;
        }
}

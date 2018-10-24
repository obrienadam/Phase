#include "CellLink.h"
#include "Cell.h"

CellLink::CellLink(const Cell &cellA, const Cell &cellB)
    :
      _cellA(&cellA),
      _cellB(&cellB),
      _rvec(cellB.centroid() - cellA.centroid())
{

}

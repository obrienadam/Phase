#include "Source.h"

namespace source
{

ScalarFiniteVolumeField div(const VectorFiniteVolumeField &field)
{
    ScalarFiniteVolumeField divF(field.grid, "divF");

    for(const Cell& cell: field.grid.cellZone("fluid"))
    {
        for(const InteriorLink& nb: cell.neighbours())
            divF(cell) += dot(field(nb.face()), nb.outwardNorm());

        for(const BoundaryLink& bd: cell.boundaries())
            divF(cell) += dot(field(bd.face()), bd.outwardNorm());
    }

    return divF;
}

}

#include "AxisymmetricSource.h"

Vector axi::src::div(const VectorFiniteVolumeField &u)
{
    Vector divU(u.grid()->localCells().size());

    for (const Cell &cell: u.cells())
    {
        Scalar tmp = 0.;
        for (const InteriorLink &nb: cell.neighbours())
            tmp += dot(u(nb.face()), nb.face().polarOutwardNorm(cell.centroid()));

        for (const BoundaryLink &bd: cell.boundaries())
            tmp += dot(u(bd.face()), bd.face().polarOutwardNorm(cell.centroid()));

        divU(u.indexMap()->local(cell, 0)) = tmp;
    }

    return divU;
}


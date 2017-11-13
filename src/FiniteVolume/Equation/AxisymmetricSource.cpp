#include "AxisymmetricSource.h"

ScalarFiniteVolumeField axi::src::div(const VectorFiniteVolumeField &u)
{
    ScalarFiniteVolumeField divU(u.gridPtr(), "", 0.);

    for(const Cell& cell: u.cells())
    {
        Scalar tmp = 0.;
        for(const InteriorLink& nb: cell.neighbours())
            tmp += dot(u(nb.face()), nb.face().polarOutwardNorm(cell.centroid()));

        for(const BoundaryLink& bd: cell.boundaries())
            tmp += dot(u(bd.face()), bd.face().polarOutwardNorm(cell.centroid()));

        divU(cell) = tmp;
    }

    return divU;
}


#include "Geometry/Tensor2D.h"

#include "Source.h"

ScalarFiniteVolumeField src::div(const VectorFiniteVolumeField& field)
{
    ScalarFiniteVolumeField divF(field.grid(), "", 0., false, false);

    for (const Cell &cell: field.cells())
    {
        Scalar div = 0.;

        for (const InteriorLink &nb: cell.neighbours())
            div += dot(field(nb.face()), nb.outwardNorm());

        for (const BoundaryLink &bd: cell.boundaries())
            div += dot(field(bd.face()), bd.outwardNorm());

        divF(cell) = div;
    }

    return divF;
}

ScalarFiniteVolumeField src::laplacian(Scalar gamma,
                                       const ScalarFiniteVolumeField &phi)
{
    ScalarFiniteVolumeField lapPhi(phi.grid(), "", 0., false, false);

    for (const Cell &cell: phi.cells())
    {
        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar coeff = gamma*dot(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();
            lapPhi(cell) += (phi(nb.cell()) - phi(cell)) * coeff;
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Scalar coeff = gamma*dot(bd.rFaceVec(), bd.outwardNorm()) / bd.rFaceVec().magSqr();
            lapPhi(cell) += (phi(bd.face()) - phi(cell)) * coeff;
        }
    }

    return lapPhi;
}

ScalarFiniteVolumeField src::laplacian(const ScalarFiniteVolumeField& gamma,
                                       const ScalarFiniteVolumeField& phi)
{
    ScalarFiniteVolumeField lapPhi(phi.grid(), "" + phi.name(), 0., false, false);

    for(const Cell& cell: phi.cells())
    {
        for (const InteriorLink& nb: cell.neighbours())
        {
            Scalar coeff = gamma(nb.face())*dot(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();
            lapPhi(cell) += (phi(nb.cell()) - phi(cell)) * coeff;
        }

        for (const BoundaryLink& bd: cell.boundaries())
        {
            Scalar coeff = gamma(bd.face()) * dot(bd.rFaceVec(), bd.outwardNorm()) / bd.rFaceVec().magSqr();
            lapPhi(cell) += (phi(bd.face()) - phi(cell)) * coeff;
        }
    }

    return lapPhi;
}

#include "Source.h"
#include "Tensor2D.h"

ScalarFiniteVolumeField src::div(const VectorFiniteVolumeField &field)
{
    ScalarFiniteVolumeField divF(field.gridPtr(), "divF", 0., false, false);

    for (const Cell &cell: field.grid().cellZone("fluid"))
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
    ScalarFiniteVolumeField lapPhi(phi.gridPtr(), "lap" + phi.name(), 0., false, false);

    for (const Cell &cell: phi.grid().cellZone("fluid"))
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
    ScalarFiniteVolumeField lapPhi(phi.gridPtr(), "lap" + phi.name(), 0., false, false);

    for(const Cell& cell: phi.grid().cellZone("fluid"))
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

VectorFiniteVolumeField src::ftc(const ScalarFiniteVolumeField& cellWeight,
                                 const ScalarFiniteVolumeField& faceWeight,
                                 const VectorFiniteVolumeField& field,
                                 const CellGroup& cells)
{
    VectorFiniteVolumeField src(field.gridPtr(), "tmp", Vector2D(0., 0.), false, false);

    for(const Cell& cell: cells)
    {
        Vector2D sumSf(0., 0.), tmp(0., 0.);

        for(const InteriorLink& nb: cell.neighbours())
        {
            Vector2D sf = nb.outwardNorm().abs();
            tmp += pointwise(field(nb.face()), sf) / faceWeight(nb.face());
            sumSf += sf;
        }

        for(const BoundaryLink& bd: cell.boundaries())
        {
            Vector2D sf = bd.outwardNorm().abs();
            tmp += pointwise(field(bd.face()), sf) / faceWeight(bd.face());
            sumSf += sf;
        }

        src(cell) = cellWeight(cell) * Vector2D(tmp.x / sumSf.x, tmp.y / sumSf.y) * cell.volume();
    }

    return src;
}

VectorFiniteVolumeField src::ftc(const ScalarFiniteVolumeField& weight,
                                 const VectorFiniteVolumeField& field,
                                 const CellGroup& cells )
{
    return ftc(weight, weight, field, cells);
}

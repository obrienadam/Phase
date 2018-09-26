#include "Geometry/Tensor2D.h"

#include "Source.h"

Vector src::div(const VectorFiniteVolumeField &field, const CellGroup &cells)
{
    Vector divU(field.grid()->localCells().size());

    for (const Cell &cell: cells)
    {
        Scalar divUc = 0.;

        for (const InteriorLink &nb: cell.neighbours())
            divUc += dot(field(nb.face()), nb.outwardNorm());

        for (const BoundaryLink &bd: cell.boundaries())
            divUc += dot(field(bd.face()), bd.outwardNorm());

        divU(field.indexMap()->local(cell, 0)) = divUc;
    }

    return divU;
}

Vector src::div(const VectorFiniteVolumeField &field)
{
    return div(field, field.cells());
}

Vector src::laplacian(Scalar gamma,
                      const ScalarFiniteVolumeField &phi)
{
    Vector lapPhi(phi.grid()->localCells().size());

    for (const Cell &cell: phi.cells())
    {
        Scalar tmp = 0.;

        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar coeff = gamma * dot(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();
            tmp += (phi(nb.cell()) - phi(cell)) * coeff;
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Scalar coeff = gamma * dot(bd.rFaceVec(), bd.outwardNorm()) / bd.rFaceVec().magSqr();
            tmp += (phi(bd.face()) - phi(cell)) * coeff;
        }

        lapPhi(phi.indexMap()->local(cell, 0)) = tmp;
    }

    return lapPhi;
}

Vector src::laplacian(const ScalarFiniteVolumeField &gamma,
                      const ScalarFiniteVolumeField &phi)
{
    Vector lapPhi(2 * phi.grid()->localCells().size());

    for (const Cell &cell: phi.cells())
    {
        Vector2D tmp = Vector2D(0., 0.);

        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar coeff = gamma(nb.face()) * dot(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();
            tmp += (phi(nb.cell()) - phi(cell)) * coeff;
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Scalar coeff = gamma(bd.face()) * dot(bd.rFaceVec(), bd.outwardNorm()) / bd.rFaceVec().magSqr();
            tmp += (phi(bd.face()) - phi(cell)) * coeff;
        }

        lapPhi(phi.indexMap()->local(cell, 0)) = tmp.x;
        lapPhi(phi.indexMap()->local(cell, 1)) = tmp.y;
    }

    return lapPhi;
}

Vector src::src(const ScalarFiniteVolumeField &field)
{
    Vector vec(field.grid()->localCells().size());

    for(const Cell& cell: field.cells())
        vec(field.indexMap()->local(cell, 0)) = field(cell) * cell.volume();

    return vec;
}

Vector src::src(const VectorFiniteVolumeField &field)
{
    Vector vec(2 * field.grid()->localCells().size());

    for(const Cell& cell: field.cells())
    {
        vec(field.indexMap()->local(cell, 0)) = field(cell).x * cell.volume();
        vec(field.indexMap()->local(cell, 1)) = field(cell).y * cell.volume();
    }

    return vec;
}

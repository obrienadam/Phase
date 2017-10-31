#include "Celeste.h"

Celeste::CelesteStencil::CelesteStencil(const Cell &cell, bool weighted)
        :
        cellPtr_(&cell)
{
    init(weighted);
}

Celeste::CelesteStencil::CelesteStencil(const Cell &cell, const ImmersedBoundary &ib, bool weighted)
        :
        cellPtr_(&cell)
{
    init(ib, weighted);
}

void Celeste::CelesteStencil::init(bool weighted)
{
    weighted_ = weighted;
    truncated_ = false;

    const Cell &cell = *cellPtr_;

    Matrix A(cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size(), 5);

    int i = 0;
    for (const InteriorLink &nb: cell.neighbours())
    {
        Vector2D r = nb.rCellVec() / (weighted_ ? nb.rCellVec().magSqr() : 1.);

        A.setRow(i++, {
                r.x * r.x / 2.,
                r.y * r.y / 2.,
                r.x * r.y,
                r.x,
                r.y
        });
    }

    for (const DiagonalCellLink &dg: cell.diagonals())
    {
        Vector2D r = dg.rCellVec() / (weighted_ ? dg.rCellVec().magSqr() : 1.);

        A.setRow(i++, {
                r.x * r.x / 2.,
                r.y * r.y / 2.,
                r.x * r.y,
                r.x,
                r.y
        });
    }

    for (const BoundaryLink &bd: cell.boundaries())
    {
        Vector2D r = bd.rFaceVec() / (weighted_ ? bd.rFaceVec().magSqr() : 1.);

        A.setRow(i++, {
                r.x * r.x / 2.,
                r.y * r.y / 2.,
                r.x * r.y,
                r.x,
                r.y
        });
    }

    pInv_ = inverse(transpose(A) * A) * transpose(A);
}

void Celeste::CelesteStencil::init(const ImmersedBoundary &ib, bool weighted)
{
    const Cell &cell = *cellPtr_;
    Matrix A(cell.neighbours().size() + cell.diagonals().size() + cell.boundaries().size(), 5);

    truncated_ = false;

    int i = 0;
    for (const InteriorLink &nb: cell.neighbours())
    {
        auto ibObj = ib.ibObj(nb.cell().centroid());

        Vector2D r = ibObj ? ibObj->intersectionLine(cell.centroid(), nb.cell().centroid()).rVec() : nb.rCellVec();

        A.setRow(i++, {
                r.x * r.x / 2.,
                r.y * r.y / 2.,
                r.x * r.y,
                r.x,
                r.y
        });

        truncated_ = ibObj ? true : truncated_;
    }

    for (const DiagonalCellLink &dg: cell.diagonals())
    {
        auto ibObj = ib.ibObj(dg.cell().centroid());

        Vector2D r = ibObj ? ibObj->intersectionLine(cell.centroid(), dg.cell().centroid()).rVec() : dg.rCellVec();

        A.setRow(i++, {
                r.x * r.x / 2.,
                r.y * r.y / 2.,
                r.x * r.y,
                r.x,
                r.y
        });

        truncated_ = ibObj ? true : truncated_;
    }

    for (const BoundaryLink &bd: cell.boundaries())
    {
        Vector2D r = bd.rFaceVec();

        A.setRow(i++, {
                r.x * r.x / 2.,
                r.y * r.y / 2.,
                r.x * r.y,
                r.x,
                r.y
        });
    }

    pInv_ = inverse(transpose(A) * A) * transpose(A);
}

Vector2D Celeste::CelesteStencil::grad(const ScalarFiniteVolumeField &phi) const
{
    Matrix b(pInv_.n(), 1);
    int i = 0;

    const Cell &cell = *cellPtr_;

    for (const InteriorLink &nb: cell.neighbours())
        b(i++, 0) = (phi(nb.cell()) - phi(cell)) / (weighted_ ? nb.rCellVec().magSqr() : 1.);

    for (const DiagonalCellLink &dg: cell.diagonals())
        b(i++, 0) = (phi(dg.cell()) - phi(cell)) / (weighted_ ? dg.rCellVec().magSqr() : 1.);

    for (const BoundaryLink &bd: cell.boundaries())
        b(i++, 0) = (phi(bd.face()) - phi(cell)) / (weighted_ ? bd.rFaceVec().magSqr() : 1.);

    b = pInv_ * b;
    return Vector2D(b(b.m() - 2, 0), b(b.m() - 1, 0));
}

Scalar Celeste::CelesteStencil::div(const VectorFiniteVolumeField &u) const
{
    Matrix b(pInv_.n(), 2);
    int i = 0;

    const Cell &cell = *cellPtr_;

    for (const InteriorLink &nb: cell.neighbours())
    {
        Vector2D du = (u(nb.cell()) - u(cell)) / (weighted_ ? nb.rCellVec().magSqr() : 1.);
        b(i, 0) = du.x;
        b(i++, 1) = du.y;
    }

    for (const DiagonalCellLink &dg: cell.diagonals())
    {
        Vector2D du = (u(dg.cell()) - u(cell)) / (weighted_ ? dg.rCellVec().magSqr() : 1.);
        b(i, 0) = du.x;
        b(i++, 1) = du.y;
    }

    for (const BoundaryLink &bd: cell.boundaries())
    {
        Vector2D du = (u(bd.face()) - u(cell)) / (weighted_ ? bd.rFaceVec().magSqr() : 1.);
        b(i, 0) = du.x;
        b(i++, 1) = du.y;
    }

    b = pInv_ * b;
    return b(b.m() - 2, 0) + b(b.m() - 1, 1);
}

Scalar Celeste::CelesteStencil::kappa(const VectorFiniteVolumeField &n,
                                      const ImmersedBoundary &ib,
                                      const Celeste &fst) const
{
    Matrix b(pInv_.n(), 2);

    const Cell &cell = *cellPtr_;

    int i = 0;
    for (const InteriorLink &nb: cell.neighbours())
    {
        auto ibObj = ib.ibObj(nb.cell().centroid());
        Vector2D dn = ((ibObj ? fst.contactLineNormal(cell, nb.cell(), *ibObj) : n(nb.cell())) - n(cell)) /
                      (weighted_ ? nb.rCellVec().magSqr() : 1.);
        b(i, 0) = dn.x;
        b(i++, 1) = dn.y;
    }

    for (const DiagonalCellLink &dg: cell.diagonals())
    {
        auto ibObj = ib.ibObj(dg.cell().centroid());
        Vector2D dn = ((ibObj ? fst.contactLineNormal(cell, dg.cell(), *ibObj) : n(dg.cell())) - n(cell)) /
                      (weighted_ ? dg.rCellVec().magSqr() : 1.);
        b(i, 0) = dn.x;
        b(i++, 1) = dn.y;
    }

    for (const BoundaryLink &bd: cell.boundaries())
    {
        Vector2D dn = (n(bd.face()) - n(cell)) / (weighted_ ? bd.rFaceVec().magSqr() : 1.);
        b(i, 0) = dn.x;
        b(i++, 1) = dn.y;
    }

    b = pInv_ * b;
    return b(b.m() - 2, 0) + b(b.m() - 1, 1);
}
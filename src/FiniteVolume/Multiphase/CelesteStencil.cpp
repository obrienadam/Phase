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

    reset();

    for (const InteriorLink &nb: cell.neighbours())
    {
        cells_.push_back(std::cref(nb.cell()));

        if (!cell.boundaries().empty())
            for (const BoundaryLink &bd: nb.cell().boundaries())
                faces_.push_back(std::cref(bd.face()));
    }

    for (const CellLink &dg: cell.diagonals())
        cells_.push_back(std::cref(dg.cell()));

    for (const BoundaryLink &bd: cell.boundaries())
        faces_.push_back(std::cref(bd.face()));

    initMatrix();
}

void Celeste::CelesteStencil::init(const ImmersedBoundary &ib, bool weighted)
{
    weighted_ = weighted;
    truncated_ = false;

    const Cell &cell = *cellPtr_;

    reset();

    for (const InteriorLink &nb: cell.neighbours())
    {
        auto ibObj = ib.ibObj(nb.cell().centroid());

        if (!ibObj)
        {
            cells_.push_back(std::cref(nb.cell()));

            if (!cell.boundaries().empty())
                for (const BoundaryLink &bd: nb.cell().boundaries())
                    faces_.push_back(std::cref(bd.face()));
        }
        else
        {
            compatPts_.push_back(
                    std::make_pair(std::cref(nb.cell()), std::weak_ptr<const ImmersedBoundaryObject>(ibObj)));
        }
    }

    for (const CellLink &dg: cell.diagonals())
    {
        auto ibObj = ib.ibObj(dg.cell().centroid());

        if (!ibObj)
            cells_.push_back(std::cref(dg.cell()));
        else
        {
            //compatPts_.push_back(std::make_pair(ibObj, dg.cell().centroid()));
            compatPts_.push_back(
                    std::make_pair(std::cref(dg.cell()), std::weak_ptr<const ImmersedBoundaryObject>(ibObj)));
        }
    }

    for (const BoundaryLink &bd: cell.boundaries())
        faces_.push_back(std::cref(bd.face()));

    initMatrix();
}

void Celeste::CelesteStencil::reset()
{
    cells_.clear();
    faces_.clear();
    compatPts_.clear();
}

Vector2D Celeste::CelesteStencil::grad(const ScalarFiniteVolumeField &phi) const
{
    Matrix b(pInv_.n(), 1);
    int i = 0;

    const Cell &cell = *cellPtr_;

    for (const Cell &kCell: cells_)
    {
        b(i++, 0) = (phi(kCell) - phi(cell)) / (weighted_ ? (kCell.centroid() - cell.centroid()).magSqr() : 1.);
    }

    for (const Face &face: faces_)
    {
        b(i++, 0) = (phi(face) - phi(cell)) / (weighted_ ? (face.centroid() - cell.centroid()).magSqr() : 1.);
    }

    b = pInv_ * b;
    return Vector2D(b(b.m() - 2, 0), b(b.m() - 1, 0));
}

Scalar Celeste::CelesteStencil::div(const VectorFiniteVolumeField &u) const
{
    Matrix b(pInv_.n(), 2);

    const Cell &cell = *cellPtr_;

    int i = 0;
    for (const Cell &kCell: cells_)
    {
        Vector2D du = (u(kCell) - u(cell)) / (weighted_ ? (kCell.centroid() - cell.centroid()).magSqr() : 1.);
        b(i, 0) = du.x;
        b(i++, 1) = du.y;
    }

    for (const Face &face: faces_)
    {
        Vector2D du = (u(face) - u(cell)) / (weighted_ ? (face.centroid() - cell.centroid()).magSqr() : 1.);
        b(i, 0) = du.x;
        b(i++, 1) = du.y;
    }

    b = pInv_ * b;
    return b(b.m() - 2, 0) + b(b.m() - 1, 1);
}

Scalar Celeste::CelesteStencil::kappa(const VectorFiniteVolumeField &n,
                                      const Celeste &fst) const
{
    return div(n);
}

Scalar Celeste::CelesteStencil::kappa(const VectorFiniteVolumeField &n,
                                      const ImmersedBoundary &ib,
                                      const Celeste &fst) const
{
    Matrix b(pInv_.n(), 2);

    const Cell &cell = *cellPtr_;

    int i = 0;
    for (const Cell &kCell: cells_)
    {
        Vector2D dn = n(kCell) - n(cell);

        if (weighted_) dn /= (kCell.centroid() - cell.centroid()).magSqr();

        b(i, 0) = dn.x;
        b(i++, 1) = dn.y;
    }

    for (const Face &face: faces_)
    {
        Vector2D dn = n(face) - n(cell);

        if (weighted_) dn /= (face.centroid() - cell.centroid()).magSqr();

        b(i, 0) = dn.x;
        b(i++, 1) = dn.y;
    }

    for (const auto &compatPt: compatPts_)
    {
        auto ibObj = compatPt.second.lock();
        Point2D pt = ibObj->intersectionLine(cell.centroid(), compatPt.first.get().centroid()).ptB();
        Vector2D dn = fst.contactLineNormal(cell, pt, *ibObj) - n(cell);

        if (weighted_) dn /= (compatPt.first.get().centroid() - cell.centroid()).magSqr();

        b(i, 0) = dn.x;
        b(i++, 1) = dn.y;
    }

    b = pInv_ * b;
    return b(b.m() - 2, 0) + b(b.m() - 1, 1);
}

//- Private methods

void Celeste::CelesteStencil::initMatrix()
{
    Matrix A(cells_.size() + compatPts_.size() + faces_.size(), 5);
    const Cell &cell = *cellPtr_;

    int i = 0;
    for (const Cell &kCell: cells_)
    {
        Vector2D r = (kCell.centroid() - cell.centroid());

        if (weighted_) r /= r.magSqr();

        A.setRow(i++, {
                r.x * r.x / 2.,
                r.y * r.y / 2.,
                r.x * r.y,
                r.x,
                r.y
        });
    }

    for (const Face &face: faces_)
    {
        Vector2D r = face.centroid() - cell.centroid();

        if (weighted_) r /= r.magSqr();

        A.setRow(i++, {
                r.x * r.x / 2.,
                r.y * r.y / 2.,
                r.x * r.y,
                r.x,
                r.y
        });
    }

    for (const auto &compatPt: compatPts_)
    {
        Vector2D r = compatPt.first.get().centroid() - cell.centroid();

        if (weighted_) r /= r.magSqr();

        A.setRow(i++, {
                r.x * r.x / 2.,
                r.y * r.y / 2.,
                r.x * r.y,
                r.x,
                r.y
        });
    }

    pInv_ = pseudoInverse(A);
}
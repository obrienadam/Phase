#include "Celeste.h"

Matrix Celeste::Stencil::b_;

Celeste::Stencil::Stencil(const Cell &cell, bool weighted)
    :
      cellPtr_(&cell)
{
    init(weighted);
}

void Celeste::Stencil::init(bool weighted)
{
    weighted_ = weighted;
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

void Celeste::Stencil::reset()
{
    cells_.clear();
    faces_.clear();
}

Vector2D Celeste::Stencil::grad(const ScalarFiniteVolumeField &phi) const
{
    b_.resize(pInv_.n(), 1);
    int i = 0;

    const Cell &cell = *cellPtr_;

    for (const Cell &kCell: cells_)
    {
        Scalar s = weighted_ ? (kCell.centroid() - cell.centroid()).magSqr() : 1.;
        b_(i++, 0) = (phi(kCell) - phi(cell)) / s;
    }

    for (const Face &face: faces_)
    {
        Scalar s = weighted_ ? (face.centroid() - cell.centroid()).magSqr() : 1.;
        b_(i++, 0) = (phi(face) - phi(cell)) / s;
    }

    b_ = pInv_ * b_;
    return Vector2D(b_(b_.m() - 2, 0), b_(b_.m() - 1, 0));
}

Scalar Celeste::Stencil::div(const VectorFiniteVolumeField &u) const
{
    b_.resize(pInv_.n(), 2);

    const Cell &cell = *cellPtr_;

    int i = 0;
    for (const Cell &kCell: cells_)
    {
        Scalar s = weighted_ ? (kCell.centroid() - cell.centroid()).magSqr() : 1.;
        Vector2D du = (u(kCell) - u(cell)) / s;
        b_(i, 0) = du.x;
        b_(i++, 1) = du.y;
    }

    for (const Face &face: faces_)
    {
        Scalar s = weighted_ ? (face.centroid() - cell.centroid()).magSqr() : 1.;
        Vector2D du = (u(face) - u(cell)) / s;
        b_(i, 0) = du.x;
        b_(i++, 1) = du.y;
    }

    b_ = pInv_ * b_;
    return b_(b_.m() - 2, 0) + b_(b_.m() - 1, 1);
}

Scalar Celeste::Stencil::axiDiv(const VectorFiniteVolumeField &u) const
{
    b_.resize(pInv_.n(), 2);

    const Cell &cell = *cellPtr_;

    int i = 0;
    for (const Cell &kCell: cells_)
    {
        Scalar s = weighted_ ? (kCell.centroid() - cell.centroid()).magSqr() : 1.;
        Vector2D du = (Vector2D(kCell.centroid().x * u(kCell).x, u(kCell).y)
                       - Vector2D(cell.centroid().x * u(cell).x, u(cell).y)) / s;
        b_(i, 0) = du.x;
        b_(i++, 1) = du.y;
    }

    for (const Face &face: faces_)
    {
        Scalar s = weighted_ ? (face.centroid() - cell.centroid()).magSqr() : 1.;
        Vector2D du = (Vector2D(face.centroid().x * u(face).x, u(face).y)
                       - Vector2D(cell.centroid().x * u(cell).x, u(cell).y)) / s;
        b_(i, 0) = du.x;
        b_(i++, 1) = du.y;
    }

    b_ = pInv_ * b_;
    return b_(b_.m() - 2, 0) / cell.centroid().x + b_(b_.m() - 1, 1);
}

Scalar Celeste::Stencil::kappa(const VectorFiniteVolumeField &n) const
{
    return div(n);
}

//- Private methods

void Celeste::Stencil::initMatrix()
{
    Matrix A(cells_.size() + faces_.size(), 5);
    const Cell &cell = *cellPtr_;

    int i = 0;
    for (const Cell &kCell: cells_)
    {
        Vector2D r = kCell.centroid() - cell.centroid();
        Scalar s = weighted_ ? r.magSqr() : 1.;

        A.setRow(i++, {
                     r.x * r.x / 2. / s,
                     r.y * r.y / 2. / s,
                     r.x * r.y / s,
                     r.x / s,
                     r.y / s
                 });
    }

    for (const Face &face: faces_)
    {
        Vector2D r = face.centroid() - cell.centroid();
        Scalar s = weighted_ ? r.magSqr() : 1.;

        A.setRow(i++, {
                     r.x * r.x / (2. * s),
                     r.y * r.y / (2. * s),
                     r.x * r.y / s,
                     r.x / s,
                     r.y / s
                 });
    }

    pInv_ = pseudoInverse(A);
}

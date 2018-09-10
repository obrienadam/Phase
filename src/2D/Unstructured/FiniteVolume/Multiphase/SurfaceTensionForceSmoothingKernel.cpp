#include "SurfaceTensionForce.h"

SurfaceTensionForce::SmoothingKernel::SmoothingKernel(const Cell &cell, Scalar eps)
    :
      cell_(cell),
      eps_(eps)
{
    kCells_ = cell_.grid().globalCells().itemsWithin(Circle(cell_.centroid(), eps));

    A_ = 0.;
    for(const Cell &kCell: kCells_)
        A_ += kernel(kCell.centroid() - cell_.centroid()) * kCell.volume();

    A_ = 1. / A_;
}

Scalar SurfaceTensionForce::SmoothingKernel::eval(const ScalarFiniteVolumeField &phi) const
{
    Scalar phiTilde = 0.;

    for(const Cell &kCell: kCells_)
        phiTilde += phi(kCell) * kernel(kCell.centroid() - cell_.centroid()) * kCell.volume();

    return A_ * phiTilde;
}

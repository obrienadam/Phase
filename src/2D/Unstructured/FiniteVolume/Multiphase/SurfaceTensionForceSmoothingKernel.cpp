#include "SurfaceTensionForce.h"

SurfaceTensionForce::SmoothingKernel::SmoothingKernel(const Cell &cell, Scalar eps, Type type)
    :
      cell_(cell),
      eps_(eps),
      type_(type)
{
    kCells_ = cell_.grid().globalCells().itemsWithin(Circle(cell_.centroid(), eps));
    setAxisymmetric(false);
}

void SurfaceTensionForce::SmoothingKernel::setAxisymmetric(bool axisymmetric)
{
    axisymmetric_ = axisymmetric;

    A_ = 0.;

    if(axisymmetric_)
        for(const Cell &kCell: kCells_)
            A_ += kernel(kCell.centroid() - cell_.centroid(), type_) * kCell.polarVolume();
    else
        for(const Cell &kCell: kCells_)
            A_ += kernel(kCell.centroid() - cell_.centroid(), type_) * kCell.volume();

    A_ = 1. / A_;
}

Scalar SurfaceTensionForce::SmoothingKernel::eval(const ScalarFiniteVolumeField &phi) const
{
    Scalar phiTilde = 0.;

    if(axisymmetric_)
        for(const Cell &kCell: kCells_)
            phiTilde += phi(kCell) * kernel(kCell.centroid() - cell_.centroid(), type_) * kCell.polarVolume();
    else
        for(const Cell &kCell: kCells_)
            phiTilde += phi(kCell) * kernel(kCell.centroid() - cell_.centroid(), type_) * kCell.volume();

    return A_ * phiTilde;
}

Scalar SurfaceTensionForce::SmoothingKernel::kernel(Vector2D dx, Type type) const
{
    switch (type)
    {
    case PESKIN:
        return kcos(dx.x) * kcos(dx.y);
    case POW_6:
    {
        Scalar r = dx.mag();
        return r < eps_ ? std::pow(std::pow(eps_, 2) - std::pow(r, 2), 3) : 0.;
    }
    case POW_8:
    {
        Scalar r = dx.mag();
        return r < eps_ ? std::pow(std::pow(eps_, 2) - std::pow(r, 2), 4) : 0.;
    }
    }
}

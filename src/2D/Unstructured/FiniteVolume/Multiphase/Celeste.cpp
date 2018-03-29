#include "Celeste.h"

Celeste::Celeste(const Input &input,
                 const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                 const std::weak_ptr<ImmersedBoundary> &ib)
        :
        SurfaceTensionForce(input, grid, ib)
{
    updateStencils();
}

void Celeste::computeFaceInterfaceForces(const ScalarFiniteVolumeField &gamma, const ScalarGradient &gradGamma)
{
    updateStencils();

    computeGradGammaTilde(gamma);
    computeInterfaceNormals();
    computeCurvature();

    auto &ft = *this;
    auto &kappa = *kappa_;

    for (const Face &face: grid_->faces())
        ft(face) = sigma_ * kappa(face) * gradGamma(face);
}

void Celeste::computeInterfaceForces(const ScalarFiniteVolumeField &gamma, const ScalarGradient &gradGamma)
{
    updateStencils();

    computeGradGammaTilde(gamma);
    computeInterfaceNormals();
    computeCurvature();

    auto &ft = *this;
    auto &kappa = *kappa_;

    ft.fill(Vector2D(0., 0.));
    for (const Cell &cell: grid_->cellZone("fluid"))
        ft(cell) = sigma_ * kappa(cell) * gradGamma(cell);

    ft.interpolateFaces();
}

//- Protected methods

void Celeste::computeGradGammaTilde(const ScalarFiniteVolumeField &gamma)
{
    smoothGammaField(gamma);

    auto &gammaTilde = *gammaTilde_;
    auto &gradGammaTilde = *gradGammaTilde_;

    gradGammaTilde.fill(Vector2D(0., 0.));
    for (const Cell &cell: gradGammaTilde.grid()->cellZone("fluid"))
        gradGammaTilde(cell) = gradGammaTildeStencils_[cell.id()].grad(gammaTilde);
}

void Celeste::computeCurvature()
{
    kappa_->fill(0.);

    auto &n = *n_;
    auto &kappa = *kappa_;
    const auto &gradGammaTilde = *gradGammaTilde_;

    for (const Cell &cell: grid_->cellZone("fluid"))
        if (gradGammaTilde(cell).magSqr() > 0.)
            kappa(cell) = kappaStencils_[cell.id()].kappa(n, *ib_.lock(), *this);

    grid_->sendMessages(kappa);

    for (const Face &face: grid_->interiorFaces())
    {
        //- According to Afkhami 2007

        if(n(face.lCell()).magSqr() != 0. && n(face.rCell()).magSqr() != 0.)
        {
            Scalar g = face.volumeWeight();
            kappa(face) = g * kappa(face.lCell()) + (1. - g) * kappa(face.rCell());
        }
        else if (n(face.lCell()).magSqr() != 0.)
            kappa(face) = kappa(face.lCell());
        else if (n(face.rCell()).magSqr() != 0.)
            kappa(face) = kappa(face.rCell());
        else
            kappa(face) = 0.;
    }
}

void Celeste::updateStencils()
{
    auto ib = ib_.lock();

    if (kappaStencils_.empty() || gradGammaTildeStencils_.empty())
    {
        kappaStencils_.resize(grid_->cells().size());
        gradGammaTildeStencils_.resize(grid_->cells().size());
        updateStencil_.assign(grid_->cells().size(), true);
    }

    if (ib)
        for (const Cell &cell: grid_->cells())
            kappaStencils_[cell.id()] = CelesteStencil(cell, *ib, false);
    else
        for (const Cell &cell: grid_->cells())
            kappaStencils_[cell.id()] = CelesteStencil(cell, false);

    for (const Cell &cell: grid_->cells())
        gradGammaTildeStencils_[cell.id()] = CelesteStencil(cell, true);
}

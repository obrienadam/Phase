#include "Celeste.h"

Celeste::Celeste(const Input &input,
                 const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                 const std::shared_ptr<CellGroup> &fluidCells)
    :
      SurfaceTensionForce(input, grid, fluidCells)
{
    computeStencils();
}

void Celeste::computeFaceInterfaceForces(const ScalarFiniteVolumeField &gamma, const ScalarGradient &gradGamma)
{
    computeGradGammaTilde(gamma);
    computeInterfaceNormals();
    computeCurvature();

    auto &fst = *fst_;
    auto &kappa = *kappa_;

    for (const Face &face: fst.grid()->faces())
        fst(face) = sigma_ * kappa(face) * gradGamma(face);
}

void Celeste::computeInterfaceForces(const ScalarFiniteVolumeField &gamma, const ScalarGradient &gradGamma)
{
    computeGradGammaTilde(gamma);
    computeInterfaceNormals();
    computeCurvature();

    auto &fst = *fst_;
    auto &kappa = *kappa_;

    fst.fill(Vector2D(0., 0.));

    for (const Cell &cell: fst.cells())
        fst(cell) = sigma_ * kappa(cell) * gradGamma(cell);
}

//- Protected methods

void Celeste::computeGradGammaTilde(const ScalarFiniteVolumeField &gamma)
{
    smoothGammaField(gamma);

    auto &gammaTilde = *gammaTilde_;
    auto &gradGammaTilde = *gradGammaTilde_;

    gradGammaTilde.fill(Vector2D(0., 0.));
    for (const Cell &cell: *fluid_)
        gradGammaTilde(cell) = gradGammaTildeStencils_[cell.id()].grad(gammaTilde);

    gradGammaTilde.sendMessages();
}

void Celeste::computeCurvature()
{
    auto &n = *n_;
    auto &kappa = *kappa_;

    auto validCurvature = [&n](const Cell &cell)
    {
        if(n(cell).magSqr() == 0.)
            return false;

        for(const CellLink &nb: cell.neighbours())
            if(n(nb.cell()).magSqr() == 0.)
                return false;

        for(const CellLink &nb: cell.diagonals())
            if(n(nb.cell()).magSqr() == 0.)
                return false;

        return true;
    };

    for (const Cell &cell: kappa.cells())
        if (validCurvature(cell))
            kappa(cell) = kappaStencils_[cell.id()].kappa(n);
        else
            kappa(cell) = 0.;

    kappa.sendMessages();

    for (const Face &face: kappa.grid()->interiorFaces())
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

    for(const Face &face: kappa.grid()->boundaryFaces())
        if(n(face.lCell()).magSqr() != 0.)
            kappa(face) = kappa(face.lCell());
}

void Celeste::computeStencils()
{
    kappaStencils_.resize(kappa_->grid()->cells().size());

    for (const Cell &cell: kappa_->grid()->cells())
        kappaStencils_[cell.id()] = Stencil(cell, false);

    gradGammaTildeStencils_.resize(gradGammaTilde_->grid()->cells().size());

    for (const Cell &cell: gradGammaTilde_->grid()->cells())
        gradGammaTildeStencils_[cell.id()] = Stencil(cell, true);
}

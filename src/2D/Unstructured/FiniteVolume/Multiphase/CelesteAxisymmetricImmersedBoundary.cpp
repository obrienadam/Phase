#include "CelesteAxisymmetricImmersedBoundary.h"

CelesteAxisymmetricImmersedBoundary::CelesteAxisymmetricImmersedBoundary(const Input &input,
                                                                         const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                                                         const std::shared_ptr<CellGroup> &fluidCells,
                                                                         const std::weak_ptr<const ImmersedBoundary> &ib)
    :
      CelesteImmersedBoundary(input, grid, fluidCells, ib)
{

}

void CelesteAxisymmetricImmersedBoundary::computeCurvature()
{
    kappa_->fill(0.);

    auto &n = *n_;
    auto &kappa = *kappa_;

    for (const Cell &cell: kappa.cells())
        if (n(cell).magSqr() > 1. - 2. * eps_ + eps_ * eps_ && n(cell).magSqr() < 1. + 2. * eps_ + eps_ * eps_)
            kappa(cell) = kappaStencils_[cell.id()].kappa(n) + 1. / cell.centroid().x;

    kappa.sendMessages();

    for (const Face &face: kappa.grid()->interiorFaces())
    {
        //- According to Afkhami 2007

        Scalar r1 = face.lCell().centroid().x;
        Scalar r2 = face.rCell().centroid().x;
        Scalar rf = face.centroid().x;

        if(rf == 0.)
            continue;

        if(kappa(face.lCell()) != 0. && kappa(face.rCell()) != 0.)
        {
            Scalar g = face.distanceWeight();
            kappa(face) = g * (kappa(face.lCell()) - 1. / r1) + (1. - g) * (kappa(face.rCell()) - 1. / r2)
                    + 1. / rf;
        }
        else if (kappa(face.lCell()) != 0.)
            kappa(face) = kappa(face.lCell()) - 1. / r1 + 1. / rf;
        else if (kappa(face.rCell()) != 0.)
            kappa(face) = kappa(face.rCell()) - 1. / r2 + 1. / rf;
        else
            kappa(face) = 0.;
    }

    for(const Face &face: kappa.grid()->boundaryFaces())
    {
        Scalar r1 = face.lCell().centroid().x;
        Scalar rf = face.centroid().x;

        if(rf == 0.)
            continue;

        if(kappa(face) != 0.)
            kappa(face) = kappa(face.lCell()) - 1. / r1 + 1. / rf;
    }
}

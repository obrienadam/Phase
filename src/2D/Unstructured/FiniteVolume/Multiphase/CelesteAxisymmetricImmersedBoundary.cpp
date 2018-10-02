#include "CelesteAxisymmetricImmersedBoundary.h"

CelesteAxisymmetricImmersedBoundary::CelesteAxisymmetricImmersedBoundary(const Input &input,
                                                                         const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                                                         const std::shared_ptr<CellGroup> &fluidCells,
                                                                         const std::weak_ptr<const ImmersedBoundary> &ib)
    :
      CelesteImmersedBoundary(input, grid, fluidCells, ib),
      kappaRZ_(std::make_shared<VectorFiniteVolumeField>(grid, "kappaRZ", 0., true, false, fluidCells))
{

}

void CelesteAxisymmetricImmersedBoundary::computeCurvature()
{
    kappa_->fill(0.);
    kappaRZ_->fill(Vector2D(0., 0.));

    auto &n = *n_;
    auto &kappa = *kappa_;

    for (const Cell &cell: kappa.cells())
        if (n(cell).magSqr() > 1. - 2. * eps_ + eps_ * eps_ && n(cell).magSqr() < 1. + 2. * eps_ + eps_ * eps_)
            kappa(cell) = kappaStencils_[cell.id()].axiDiv(n);

    kappa.sendMessages();

    for (const Face &face: kappa.grid()->interiorFaces())
    {
        //- According to Afkhami 2007

        if(kappa(face.lCell()) != 0. && kappa(face.rCell()) != 0.)
        {
            Scalar g = face.distanceWeight();
            kappa(face) = g * (kappa(face.lCell())) + (1. - g) * (kappa(face.rCell()));
        }
        else if (kappa(face.lCell()) != 0.)
            kappa(face) = kappa(face.lCell());
        else if (kappa(face.rCell()) != 0.)
            kappa(face) = kappa(face.rCell());
        else
            kappa(face) = 0.;
    }

    for(const Face &face: kappa.grid()->boundaryFaces())
        if(kappa(face.lCell()) != 0.)
            kappa(face) = kappa(face.lCell());
}

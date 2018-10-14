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
    auto &n = *n_;
    auto &kappa = *kappa_;

    auto validCurvature = [&n](const Cell &c)
    {
        if(n(c).magSqr() == 0.)
            return false;

        for(const CellLink &nb: c.neighbours())
            if(n(nb.cell()).magSqr() == 0.)
                return false;

        for(const CellLink &nb: c.diagonals())
            if(n(nb.cell()).magSqr() == 0.)
                return false;

        return true;
    };

    for (const Cell &cell: kappa.cells())
        if (validCurvature(cell))
            kappa(cell) = kappaStencils_[cell.id()].axiDiv(n);
        else
            kappa(cell) = 0.;

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

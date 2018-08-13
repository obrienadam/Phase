#include "CelesteImmersedBoundary.h"
#include "Geometry/Intersection.h"

CelesteImmersedBoundary::CelesteImmersedBoundary(const Input &input,
                                                 const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                                 const std::shared_ptr<CellGroup> &fluidCells,
                                                 const std::weak_ptr<const ImmersedBoundary> &ib)
    :
      Celeste(input, grid, fluidCells),
      ib_(ib)
{
    for(const auto& ibObj: *ib_.lock())
    {
        ibContactAngles_[ibObj->name()] = input.boundaryInput().get<Scalar>(
                    "ImmersedBoundaries." + ibObj->name() + ".gamma.contactAngle",
                    90.) * M_PI / 180.;
    }
}

Scalar CelesteImmersedBoundary::theta(const ImmersedBoundaryObject &ibObj) const
{
    auto it = ibContactAngles_.find(ibObj.name());
    return it != ibContactAngles_.end() ? it->second : M_PI_2;
}

void CelesteImmersedBoundary::computeContactLineExtension(ScalarFiniteVolumeField &gamma) const
{
    for(const std::shared_ptr<ImmersedBoundaryObject> &ibObj: *ib_.lock())
        for(const Cell& cell: ibObj->cells())
            if(ibObj->isInIb(cell.centroid()))
            {
                Scalar distSqr = (ibObj->nearestIntersect(cell.centroid()) - cell.centroid()).magSqr();

                if(distSqr <= kernelWidth_ * kernelWidth_)
                {
                    ContactLineStencil st(cell,
                                          *ibObj,
                                          ibContactAngles_.find(ibObj->name())->second,
                                          gamma);

                    gamma(cell) = st.gamma();
                }
            }
}

void CelesteImmersedBoundary::contactLineBcs(FiniteVolumeEquation<Scalar> &gammaEqn) const
{
    const ScalarFiniteVolumeField &gamma = gammaEqn.field();

    for(const std::shared_ptr<ImmersedBoundaryObject> &ibObj: *ib_.lock())
        for(const Cell& cell: ibObj->cells())
        {
            Scalar theta = ibContactAngles_.find(ibObj->name())->second;

            if(ibObj->isInIb(cell.centroid()))
            {
                gammaEqn.remove(cell);

                Scalar distSqr = (ibObj->nearestIntersect(cell.centroid()) - cell.centroid()).magSqr();

                if(distSqr <= kernelWidth_ * kernelWidth_)
                {
                    ContactLineStencil st(cell,
                                          *ibObj,
                                          theta,
                                          gamma);

                    Scalar alpha = st.link().alpha(st.cl()[2]);

                    gammaEqn.add(cell, cell, -1.);
                    gammaEqn.add(cell, st.link().self(), alpha);
                    gammaEqn.add(cell, st.link().cell(), 1. - alpha);
                }
                else
                {
                    gammaEqn.remove(cell);
                    gammaEqn.add(cell, cell, -1.);
                    gammaEqn.addSource(cell, 0.);
                }
            }
        }
}

void CelesteImmersedBoundary::computeInterfaceNormals()
{
    const VectorFiniteVolumeField &gradGammaTilde = *gradGammaTilde_;
    VectorFiniteVolumeField &n = *n_;

    for (const Cell &cell: n.grid()->cells())
        n(cell) = gradGammaTilde(cell).magSqr() >= eps_ * eps_ ? -gradGammaTilde(cell).unitVec() : Vector2D(0., 0.);

    //- Override the ib cells
    for(const auto &ibObj: *ib_.lock())
        for(const Cell& cell: ibObj->cells())
            if(ibObj->isInIb(cell.centroid()))
            {
                Scalar distSqr = (ibObj->nearestIntersect(cell.centroid()) - cell.centroid()).magSqr();

                if(distSqr <= kernelWidth_ * kernelWidth_)
                {
                    ContactLineStencil st(cell,
                                          *ibObj,
                                          ibContactAngles_.find(ibObj->name())->second,
                                          *gammaTilde_);

                    n(cell) = st.ncl();
                }
            }

    //- Boundary faces set from contact line orientation
    for (const FaceGroup &patch: n.grid()->patches())
    {
        Scalar theta = SurfaceTensionForce::theta(patch);

        for (const Face &face: patch)
        {
            if (n(face.lCell()).magSqr() == 0.)
            {
                n(face) = Vector2D(0., 0.);
                continue;
            }

            Vector2D ns = -face.outwardNorm(face.lCell().centroid()).unitVec();
            Vector2D ts = (n(face.lCell()) - dot(n(face.lCell()), ns) * ns).unitVec();

            n(face) = ns * std::cos(theta) + ts * std::sin(theta);
        }
    }

    grid_->sendMessages(n);
}

void CelesteImmersedBoundary::computeCurvature()
{
    kappa_->fill(0.);

    auto &n = *n_;
    auto &kappa = *kappa_;

    for (const Cell &cell: kappa.cells())
    {
        bool isInterface = n(cell).magSqr() > 0. && !ib_.lock()->ibObj(cell.centroid());

        if(isInterface)
            for(const CellLink &nb: cell.cellLinks())
                if(n(nb.cell()).magSqr() == 0.)
                {
                    isInterface = false;
                    break;
                }

        if (isInterface)
            kappa(cell) = kappaStencils_[cell.id()].kappa(n);
    }

    kappa.grid()->sendMessages(kappa);

    for (const Face &face: kappa.grid()->interiorFaces())
    {
        //- According to Afkhami 2007

        if(kappa(face.lCell()) != 0. && kappa(face.rCell()) != 0.)
        {
            Scalar g = face.volumeWeight();
            kappa(face) = g * kappa(face.lCell()) + (1. - g) * kappa(face.rCell());
        }
        else if (kappa(face.lCell()) != 0.)
            kappa(face) = kappa(face.lCell());
        else if (kappa(face.rCell()) != 0.)
            kappa(face) = kappa(face.rCell());
        else
            kappa(face) = 0.;
    }

    for(const Face &face: kappa.grid()->boundaryFaces())
        if(n(face.lCell()).magSqr() != 0.)
            kappa(face) = kappa(face.lCell());
}

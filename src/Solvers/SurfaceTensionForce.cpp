#include "SurfaceTensionForce.h"

SurfaceTensionForce::SurfaceTensionForce(const Input &input,
                                         const ImmersedBoundary &ib,
                                         const ScalarFiniteVolumeField& gamma,
                                         const ScalarFiniteVolumeField& rho,
                                         const ScalarFiniteVolumeField& mu,
                                         const VectorFiniteVolumeField& u,
                                         const ScalarGradient& gradGamma)
        :
        VectorFiniteVolumeField(gamma.gridPtr(), "ft", Vector2D(0., 0.), true, false),
        ib_(ib),
        gamma_(gamma),
        rho_(rho),
        mu_(mu),
        u_(u),
        gradGamma_(gradGamma),
        kappa_(std::make_shared<ScalarFiniteVolumeField>(grid_, "kappa")),
        gammaTilde_(std::make_shared<ScalarFiniteVolumeField>(grid_, "gammaTilde")),
        gradGammaTilde_(std::make_shared<ScalarGradient>(*gammaTilde_)),
        n_(std::make_shared<VectorFiniteVolumeField>(grid_, "n"))
{
    sigma_ = input.caseInput().get<Scalar>("Properties.sigma");
    kernelWidth_ = input.caseInput().get<Scalar>("Solver.smoothingKernelRadius");

    //- Must figure out which patches to enforce the contact angle on
    for (const auto &input: input.boundaryInput().get_child("Boundaries." + gamma.name()))
    {
        patchContactAngles_.insert(std::make_pair(
                grid_->patch(input.first).id(),
                input.second.get<Scalar>("contactAngle", 90) * M_PI / 180.
        ));
    }

    //- Set contact angle for ibs
    for (const auto &input: input.boundaryInput().get_child("ImmersedBoundaries"))
    {
        ibContactAngles_.insert(std::make_pair(
                ib_.ibObj(input.first).id(),
                input.second.get<Scalar>(gamma.name() + ".contactAngle", 90) * M_PI / 180.
        ));
    }
}

void SurfaceTensionForce::computeInterfaceNormals()
{
    n_->compute([this](const Cell& cell){
        return gradGamma_(cell).magSqr()  > 0. ? -gradGamma_(cell).unitVec() : Vector2D(0., 0.);
    });

    n_->compute([this](const Face& face){
        return gradGamma_(face).magSqr() > 0. ? -gradGamma_(face).unitVec() : Vector2D(0., 0.);
    });
}

Vector2D SurfaceTensionForce::contactLineNormal(const Face& face)
{
    Vector2D wn = face.outwardNorm(face.lCell().centroid()).unitVec();
    Scalar theta = this->theta(grid_->patch(face));
    return dot(-gradGamma_(face.lCell()), wn) > 0. ? wn.rotate(M_PI / 2. - theta) : wn.rotate(M_PI / 2. + theta);
}

Vector2D SurfaceTensionForce::contactLineNormal(const Cell& lCell,
                                                const Cell& rCell,
                                                const ImmersedBoundaryObject& ibObj)
{
    LineSegment2D ln = ibObj.intersectionLine(lCell.centroid(), rCell.centroid());
}
#include "SurfaceTensionForce.h"

SurfaceTensionForce::SurfaceTensionForce(const Input &input,
                                         const ScalarFiniteVolumeField& gamma,
                                         const ScalarFiniteVolumeField& rho,
                                         const ScalarFiniteVolumeField& mu,
                                         const VectorFiniteVolumeField& u,
                                         VectorFiniteVolumeField& gradGamma)
        :
    VectorFiniteVolumeField(gamma.gridPtr(), "ft", Vector2D(0., 0.), true, false),
    gamma_(gamma),
    rho_(rho),
    mu_(mu),
    u_(u),
    gradGamma_(gradGamma),
    kappa_(std::make_shared<ScalarFiniteVolumeField>(grid_, "kappa")),
    n_(std::make_shared<VectorFiniteVolumeField>(grid_, "n"))
{
    sigma_ = input.caseInput().get<Scalar>("Properties.sigma");
    thetaAdv_ = input.caseInput().get<Scalar>("Properties.advancingContactAngle") * M_PI / 180.;
    thetaRec_ = input.caseInput().get<Scalar>("Properties.recedingContactAngle") * M_PI / 180.;

    std::string contactAngleType = input.caseInput().get<std::string>("Solver.contactAngleType", "static");

    //- Must figure out which patches to enforce the contact angle on
    for (const auto &boundaryInput: input.boundaryInput().get_child("Boundaries.gamma"))
    {
        std::string status = boundaryInput.second.get<std::string>("contactAngle", "off");

        if (status == "on")
        {
            contactAnglePatches_.push_back(gamma_.grid().patch(boundaryInput.first));
            printf("Contact angles will be applied to patch \"%s\".\n", boundaryInput.first.c_str());
        }
        else if (status == "off")
            continue;
        else
            throw Exception("SurfaceTensionForce", "SurfaceTensionForce",
                            "invalid contact angle status \"" + status + "\".");
    }
}

Scalar SurfaceTensionForce::theta(const Cell& cell)
{
    return thetaAdv_;
}

Vector2D SurfaceTensionForce::computeContactLineNormal(const Vector2D &gradGamma,
                                                       const Vector2D &wallNormal,
                                                       const Vector2D &vel,
                                                       Scalar theta) const
{
    Vector2D nt = wallNormal.tangentVec();
    return dot(-gradGamma, nt) > 0. ? nt.rotate(M_PI / 2. - theta).unitVec() : nt.rotate(M_PI / 2. + theta).unitVec();
}

Vector2D SurfaceTensionForce::computeContactLineNormal(const Vector2D &gradGamma,
                                                       const Vector2D &wallNormal,
                                                       const Vector2D &vel) const
{
    Vector2D nt = wallNormal.tangentVec();
    Scalar theta = dot(-gradGamma, vel) >= 0. ? thetaAdv_ : thetaRec_;
    return dot(-gradGamma, nt) > 0. ? nt.rotate(M_PI / 2. - theta).unitVec() : nt.rotate(M_PI / 2. + theta).unitVec();
}

void SurfaceTensionForce::computeGhostCellVals(const ImmersedBoundaryObject &ibObj, Scalar theta)
{
    throw Exception("SurfaceTensionForce", "computeGhostCellVals", "deprecated.");
}

bool SurfaceTensionForce::isContactLinePatch(const Patch &patch) const
{
    for (const Patch &clPatch: contactAnglePatches_)
    {
        if (clPatch.id() == patch.id())
            return true;
    }

    return false;
}

Scalar SurfaceTensionForce::ibTheta(const ImmersedBoundaryObject &ibObj) const
{
    return ibContactAngles_.find(ibObj.id())->second;
}

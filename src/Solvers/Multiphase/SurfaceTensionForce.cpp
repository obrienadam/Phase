#include "SurfaceTensionForce.h"

SurfaceTensionForce::SurfaceTensionForce(const Input &input, const ScalarFiniteVolumeField &gamma)
    :
      gamma_(gamma),
      n_(gamma.grid, "interfaceNormals"),
      kappa_(gamma.grid, "interfaceCurvature")
{
    sigma_ = input.caseInput().get<Scalar>("Properties.sigma");
    thetaAdv_ = input.caseInput().get<Scalar>("Properties.advancingContactAngle", 90)*M_PI/180.;
    thetaRec_ = input.caseInput().get<Scalar>("Properties.recedingContactAngle", 90)*M_PI/180.;

    //- Must figure out which patches to enforce the contact angle on
    for(const auto &boundaryInput: input.boundaryInput().get_child("Boundaries.gamma"))
    {
        std::string status = boundaryInput.second.get<std::string>("contactAngle", "off");
/*
        if(status == "on")
            contactAnglePatches_.push_back(std::cref(gamma_.grid.patches().find(boundaryInput.first)->second));
        else if(status == "off")
            continue;
        else
            throw Exception("SurfaceTensionForce", "SurfaceTensionForce", "invalid contact angle status \"" + status + "\".");*/
    }
}


Vector2D SurfaceTensionForce::computeContactLineNormal(const Vector2D& gradGamma, const Vector2D& wallNormal) const
{
    Vector2D nt = wallNormal.rotate(M_PI/2.);

    return dot(-gradGamma, nt) > 0. ? nt.rotate(M_PI/2. -thetaAdv_).unitVec() : nt.rotate(M_PI/2. + thetaAdv_).unitVec();
}

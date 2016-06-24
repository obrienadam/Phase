#include "SurfaceTensionForce.h"

SurfaceTensionForce::SurfaceTensionForce(const Input &input, const ScalarFiniteVolumeField &gamma, const VectorFiniteVolumeField &u, std::map<std::string, VectorFiniteVolumeField> &fields)
    :
      gamma_(gamma),
      u_(u),
      n_((fields.insert(std::make_pair(std::string("n"), VectorFiniteVolumeField(gamma.grid, "n"))).first)->second),
      kappa_(gamma.grid, "interfaceCurvature")
{
    sigma_ = input.caseInput().get<Scalar>("Properties.sigma");
    thetaAdv_ = input.caseInput().get<Scalar>("Properties.advancingContactAngle", 90)*M_PI/180.;
    thetaRec_ = input.caseInput().get<Scalar>("Properties.recedingContactAngle", 90)*M_PI/180.;

    //- Must figure out which patches to enforce the contact angle on
    for(const auto &boundaryInput: input.boundaryInput().get_child("Boundaries.gamma"))
    {
        std::string status = boundaryInput.second.get<std::string>("contactAngle", "off");

        if(status == "on")
            contactAnglePatches_.push_back(std::cref(gamma_.grid.patches().find(boundaryInput.first)->second));
        else if(status == "off")
            continue;
        else
            throw Exception("SurfaceTensionForce", "SurfaceTensionForce", "invalid contact angle status \"" + status + "\".");
    }
}


Vector2D SurfaceTensionForce::computeContactLineNormal(const Vector2D& gradGamma, const Vector2D& wallNormal, const Vector2D& vel) const
{
    const Vector2D nt = wallNormal.tangentVec();
    const Scalar theta = dot(-gradGamma, vel) >= 0. ? thetaAdv_ : thetaRec_;
    return dot(-gradGamma, nt) > 0. ? nt.rotate(M_PI/2. - theta).unitVec() : nt.rotate(M_PI/2. + theta).unitVec();
}

#include "SurfaceTensionForce.h"

SurfaceTensionForce::SurfaceTensionForce(const Input &input,
                                         Solver &solver)
    :
      gamma_(solver.scalarFields().find("gamma")->second),
      u_(solver.vectorFields().find("u")->second),
      n_(solver.addVectorField("n")),
      kappa_(solver.addScalarField("kappa")),
      gradGamma_(solver.vectorFields().find("gradGamma")->second),
      solver_(solver)
{
    sigma_ = input.caseInput().get<Scalar>("Properties.sigma");
    thetaAdv_ = input.caseInput().get<Scalar>("Properties.advancingContactAngle")*M_PI/180.;
    thetaRec_ = input.caseInput().get<Scalar>("Properties.recedingContactAngle")*M_PI/180.;

    //- Must figure out which patches to enforce the contact angle on
    for(const auto &boundaryInput: input.boundaryInput().get_child("Boundaries.gamma"))
    {
        std::string status = boundaryInput.second.get<std::string>("contactAngle", "off");

        if(status == "on")
        {
            contactAnglePatches_.push_back(std::cref(gamma_.grid.patches().find(boundaryInput.first)->second));
            printf("Contact angles will be applied to patch \"%s\".\n", boundaryInput.first.c_str());
        }
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

bool SurfaceTensionForce::isContactLinePatch(const Patch &patch) const
{
    for(const Patch &clPatch: contactAnglePatches_)
    {
        if(clPatch.id() == patch.id())
            return true;
    }

    return false;
}

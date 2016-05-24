#include "SurfaceTensionForce.h"

SurfaceTensionForce::SurfaceTensionForce(const Input &input, const ScalarFiniteVolumeField &gamma)
    :
      gamma_(gamma)
{
    sigma_ = input.caseInput().get<Scalar>("Properties.sigma");
    thetaAdv_ = input.caseInput().get<Scalar>("Properties.advancingContactAngle", 90)*M_PI/180.;
    thetaRec_ = input.caseInput().get<Scalar>("Properties.recedingContactAngle", 90)*M_PI/180.;
}

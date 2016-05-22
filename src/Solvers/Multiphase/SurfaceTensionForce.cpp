#include "SurfaceTensionForce.h"

SurfaceTensionForce::SurfaceTensionForce(const Input &input, const ScalarFiniteVolumeField &gamma)
    :
      gamma_(gamma)
{
    sigma_ = input.caseInput().get<Scalar>("Properties.sigma");
}

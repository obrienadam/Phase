#include "HeightFunction.h"

HeightFunction::HeightFunction(const Input &input, const ScalarFiniteVolumeField &gamma, const VectorFiniteVolumeField& u)
    :
      SurfaceTensionForce(input, gamma, u)
{

}

VectorFiniteVolumeField HeightFunction::compute()
{

}

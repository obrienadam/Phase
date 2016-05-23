#include "HeightFunction.h"

HeightFunction::HeightFunction(const Input &input, const ScalarFiniteVolumeField &gamma)
    :
      SurfaceTensionForce(input, gamma)
{

}

VectorFiniteVolumeField HeightFunction::compute()
{

}

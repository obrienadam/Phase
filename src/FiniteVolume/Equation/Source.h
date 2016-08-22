#ifndef SOURCE_H
#define SOURCE_H

#include "VectorFiniteVolumeField.h"

namespace fv
{

VectorFiniteVolumeField source(VectorFiniteVolumeField field);
VectorFiniteVolumeField gravity(const ScalarFiniteVolumeField &rho, const Vector2D& g);
ScalarFiniteVolumeField hydroStaticPressureBoundaries(const ScalarFiniteVolumeField& p, const ScalarFiniteVolumeField& rho, const Vector2D& g);

}

#endif

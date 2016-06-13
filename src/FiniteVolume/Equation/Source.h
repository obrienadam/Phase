#ifndef SOURCE_H
#define SOURCE_H

#include "VectorFiniteVolumeField.h"

namespace fv
{

VectorFiniteVolumeField grad(const ScalarFiniteVolumeField& field);
VectorFiniteVolumeField source(VectorFiniteVolumeField field);
VectorFiniteVolumeField gravity(const ScalarFiniteVolumeField &rho, const Vector2D& g);

}

#endif

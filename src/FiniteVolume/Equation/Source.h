#ifndef SOURCE_H
#define SOURCE_H

#include "VectorFiniteVolumeField.h"

namespace fv
{

VectorFiniteVolumeField grad(const ScalarFiniteVolumeField& field);
VectorFiniteVolumeField source(VectorFiniteVolumeField field);

}

#endif

#ifndef SOURCE_H
#define SOURCE_H

#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

namespace source
{

    ScalarFiniteVolumeField div(const VectorFiniteVolumeField &field);

    VectorFiniteVolumeField grad(const ScalarFiniteVolumeField &field);

}

#endif

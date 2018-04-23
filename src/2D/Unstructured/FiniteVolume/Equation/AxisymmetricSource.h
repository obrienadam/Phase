#ifndef PHASE_AXISYMMETRIC_SOURCE_H
#define PHASE_AXISYMMETRIC_SOURCE_H

#include "Math/Vector.h"

#include "FiniteVolume/Field/VectorFiniteVolumeField.h"

namespace axi
{
    namespace src
    {
        Vector div(const VectorFiniteVolumeField &u);
    }
}

#endif

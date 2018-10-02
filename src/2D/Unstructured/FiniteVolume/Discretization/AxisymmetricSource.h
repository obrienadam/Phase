#ifndef PHASE_AXISYMMETRIC_SOURCE_H
#define PHASE_AXISYMMETRIC_SOURCE_H

#include "Math/Vector.h"

#include "FiniteVolume/Field/VectorFiniteVolumeField.h"

namespace axi
{
    namespace src
    {
        Vector src(const ScalarFiniteVolumeField &phi);

        Vector src(const VectorFiniteVolumeField &u);

        Vector div(const VectorFiniteVolumeField &u);
    }
}

#endif

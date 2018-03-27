#ifndef AXISYMMETRIC_SOURCE_H
#define AXISYMMETRIC_SOURCE_H

#include "FiniteVolume/Field/ScalarFiniteVolumeField.h"
#include "FiniteVolume/Field/VectorFiniteVolumeField.h"

namespace axi
{
    namespace src
    {
        ScalarFiniteVolumeField div(const VectorFiniteVolumeField &u);
    }
}

#endif

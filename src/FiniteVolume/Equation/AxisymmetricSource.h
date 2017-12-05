#ifndef AXISYMMETRIC_SOURCE_H
#define AXISYMMETRIC_SOURCE_H

#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

namespace axi::src
{
    ScalarFiniteVolumeField div(const VectorFiniteVolumeField &u);
}

#endif

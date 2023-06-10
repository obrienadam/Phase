#ifndef PHASE_SOURCE_H
#define PHASE_SOURCE_H

#include "Math/Vector.h"

#include "FiniteVolume/Field/ScalarFiniteVolumeField.h"
#include "FiniteVolume/Field/VectorFiniteVolumeField.h"

namespace src {
Vector div(const VectorFiniteVolumeField &field, const CellGroup &cells);

Vector div(const VectorFiniteVolumeField &field);

Vector laplacian(Scalar gamma, const ScalarFiniteVolumeField &phi);

Vector laplacian(const ScalarFiniteVolumeField &gamma,
                 const ScalarFiniteVolumeField &phi);

Vector src(const ScalarFiniteVolumeField &field);

Vector src(const VectorFiniteVolumeField &field);
} // namespace src

#endif

#ifndef VECTOR_FINITE_VOLUME_FIELD_H
#define VECTOR_FINITE_VOLUME_FIELD_H

#include "FiniteVolumeField.h"
#include "Vector2D.h"
#include "ScalarFiniteVolumeField.h"

typedef FiniteVolumeField<Vector2D> VectorFiniteVolumeField;

VectorFiniteVolumeField grad(const ScalarFiniteVolumeField& scalarField);

VectorFiniteVolumeField operator*(const ScalarFiniteVolumeField& lhs, VectorFiniteVolumeField rhs);
VectorFiniteVolumeField operator*(VectorFiniteVolumeField lhs, const ScalarFiniteVolumeField& rhs);
VectorFiniteVolumeField operator*(const ScalarFiniteVolumeField& lhs, const Vector2D& rhs);

#endif

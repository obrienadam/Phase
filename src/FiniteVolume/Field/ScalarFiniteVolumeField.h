#ifndef SCALAR_FINITE_VOLUME_FIELD_H
#define SCALAR_FINITE_VOLUME_FIELD_H

#include "FiniteVolumeField.h"

typedef FiniteVolumeField<Scalar> ScalarFiniteVolumeField;

//- Specializations
template<>
ScalarFiniteVolumeField &ScalarFiniteVolumeField::operator=(const Vector &rhs);

template<>
void ScalarFiniteVolumeField::setBoundaryRefValues(const Input &input);

//- External
ScalarFiniteVolumeField operator*(const ScalarFiniteVolumeField &lhs, ScalarFiniteVolumeField rhs);

ScalarFiniteVolumeField operator/(ScalarFiniteVolumeField lhs, const ScalarFiniteVolumeField &rhs);

ScalarFiniteVolumeField operator/(Scalar lhs, ScalarFiniteVolumeField rhs);

#endif

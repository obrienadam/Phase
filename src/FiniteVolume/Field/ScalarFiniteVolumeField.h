#ifndef SCALAR_FINITE_VOLUME_FIELD_H
#define SCALAR_FINITE_VOLUME_FIELD_H

#include "FiniteVolumeField.h"

class ScalarFiniteVolumeField : public FiniteVolumeField<Scalar>
{
public:
    using FiniteVolumeField<Scalar>::FiniteVolumeField;
};

#endif

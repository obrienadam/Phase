#ifndef SCALAR_FINITE_VOLUME_FIELD_H
#define SCALAR_FINITE_VOLUME_FIELD_H

#include "FiniteVolumeField.h"
#include "SparseVector.h"

class ScalarFiniteVolumeField : public FiniteVolumeField<Scalar>
{
public:
    using FiniteVolumeField<Scalar>::FiniteVolumeField;

    ScalarFiniteVolumeField& operator =(const SparseVector& rhs);
};

#endif

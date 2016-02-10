#ifndef VECTOR_FINITE_VOLUME_FIELD_H
#define VECTOR_FINITE_VOLUME_FIELD_H

#include "FiniteVolumeField.h"

class VectorFiniteVolumeField : public FiniteVolumeField<Vector2D>
{
public:
    using FiniteVolumeField<Vector2D>::FiniteVolumeField;
};

#endif

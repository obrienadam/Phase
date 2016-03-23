#ifndef VECTOR_FINITE_VOLUME_FIELD_H
#define VECTOR_FINITE_VOLUME_FIELD_H

#include "Field.h"
#include "Vector2D.h"

class FiniteVolumeGrid2D;

class VectorFiniteVolumeField : public Field<Vector2D>
{
public:

    VectorFiniteVolumeField(const FiniteVolumeGrid2D& grid, std::string name);

protected:

    std::vector<Vector2D> faces_;

};

#endif

#ifndef PHASE_COLLISION_MODEL_H
#define PHASE_COLLISION_MODEL_H

#include "ImmersedBoundaryObject.h"

class CollisionModel
{
public:

    CollisionModel(Scalar eps, Scalar range = 0.);

    virtual Vector2D force(const ImmersedBoundaryObject& ibObjP, const ImmersedBoundaryObject& ibObjQ) const;

    virtual Vector2D force(const ImmersedBoundaryObject& ibObj, const FiniteVolumeGrid2D& grid) const;

private:

    Scalar eps_, range_;
};


#endif

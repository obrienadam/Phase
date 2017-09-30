#ifndef COLLISION_MODEL_H
#define COLLISION_MODEL_H

#include "ImmersedBoundaryObject.h"

class CollisionModel
{
public:

    CollisionModel(Scalar eps, Scalar range = 0.);

    virtual Vector2D force(const ImmersedBoundaryObject& ibObjP, const ImmersedBoundaryObject& ibObjQ) const;

private:

    Scalar eps_, range_;
};


#endif

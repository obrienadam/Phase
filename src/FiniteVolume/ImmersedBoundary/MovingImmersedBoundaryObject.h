#ifndef MOVING_IMMERSED_BOUNDARY_OBJECT_H
#define MOVING_IMMERSED_BOUNDARY_OBJECT_H

#include "ImmersedBoundaryObject.h"

class MovingImmersedBoundaryObject : public ImmersedBoundaryObject
{
public:
    MovingImmersedBoundaryObject(const ImmersedBoundaryObject& ibObj);

    virtual void update(Scalar timeStep) = 0;

protected:

    void updateCells();
};

#include "TranslatingImmersedBoundaryObject.h"

#endif

#ifndef TRANSLATING_IMMERSED_BOUNDARY_OBJECT
#define TRANSLATING_IMMERSED_BOUNDARY_OBJECT

#include "MovingImmersedBoundaryObject.h"

class TranslatingImmersedBoundaryObject: public MovingImmersedBoundaryObject
{
public:
    void update(Scalar timeStep);

private:

    Vector2D velocity_;
};

#endif

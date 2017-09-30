#ifndef TRANSLATING_MOTION_H
#define TRANSLATING_MOTION_H

#include "Motion.h"
#include "ImmersedBoundaryObject.h"

class TranslatingMotion: public Motion
{
public:

    TranslatingMotion(std::weak_ptr<ImmersedBoundaryObject> ibObj,
                      const Vector2D& vel,
                      const Vector2D& acc = Vector2D(0., 0.));

    void update(Scalar timeStep);
};

#endif
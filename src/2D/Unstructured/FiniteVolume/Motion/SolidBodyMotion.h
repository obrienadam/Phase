#ifndef SOLID_BODY_MOTION_H
#define SOLID_BODY_MOTION_H

#include "FiniteVolume/ImmersedBoundary/ImmersedBoundaryObject.h"

#include "Motion.h"

class SolidBodyMotion: public Motion
{
public:

    SolidBodyMotion(std::weak_ptr<ImmersedBoundaryObject> ibObj, const Vector2D& v0 = Vector2D(0., 0.));

    void update(Scalar timeStep);

private:

    Vector2D force_;
    Scalar torque_;

    std::weak_ptr<ImmersedBoundaryObject> ibObj_;
};

#endif

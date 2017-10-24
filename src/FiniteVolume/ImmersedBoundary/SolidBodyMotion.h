#ifndef SOLID_BODY_MOTION_H
#define SOLID_BODY_MOTION_H

#include "Motion.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

class SolidBodyMotion: public Motion
{
public:

    SolidBodyMotion(std::weak_ptr<ImmersedBoundaryObject> ibObj);

    void update(Scalar timeStep);

private:

    Vector2D force_;
    Scalar torque_;

};

#endif
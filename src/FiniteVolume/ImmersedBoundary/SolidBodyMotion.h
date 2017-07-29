#ifndef SOLID_BODY_MOTION_H
#define SOLID_BODY_MOTION_H

#include "Motion.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

class SolidBodyMotion: public Motion
{
public:

    SolidBodyMotion(const Vector2D& pos,
                    Scalar materialDensity);

    void update(ImmersedBoundaryObject& ibObj, Scalar timeStep);

private:

    Scalar materialDensity_; // object properties

    Vector2D force_, force0_, g_;
    Scalar torque_, torque0_;
};

#endif
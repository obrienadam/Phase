#include "SolidBodyMotion.h"
#include "ImmersedBoundaryObject.h"

SolidBodyMotion::SolidBodyMotion(const Vector2D& pos,
                                 Scalar materialDensity)
        :
        Motion(pos),
        materialDensity_(materialDensity)
{

}

void SolidBodyMotion::update(ImmersedBoundaryObject &ibObj, Scalar timeStep)
{
    force0_ = force_;
    force_ = ibObj.force();
    Scalar mass = ibObj.shape().area() * materialDensity_;

    vel_ += g_*timeStep + timeStep/(2*mass) * (force0_ + force_);

    torque0_ = torque_;
    torque_ = ibObj.torque();

    Scalar I = ibObj.shape().momentOfInertia() * materialDensity_;

    omega_ += timeStep/(2*I) * (torque0_ + torque_);
}
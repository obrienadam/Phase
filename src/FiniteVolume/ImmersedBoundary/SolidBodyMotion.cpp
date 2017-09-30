#include "SolidBodyMotion.h"
#include "ImmersedBoundaryObject.h"

SolidBodyMotion::SolidBodyMotion(std::weak_ptr<ImmersedBoundaryObject> ibObj)
        :
        Motion(ibObj)
{
    if(ibObj.lock()->mass() <= 0.)
        throw Exception("SolidBodyMotion",
                        "SolidBodyMotion",
                        "must specify a non-zero density for immersed boundary object \"" + ibObj.lock()->name() + "\".");
}

void SolidBodyMotion::update(Scalar timeStep)
{
    auto ibObj = ibObj_.lock();

    Vector2D f0 = force_;
    force_ = ibObj->force();

    Vector2D v0 = vel_;
    vel_ += timeStep / (2. * ibObj->mass()) * (force_ + f0);
    pos_ += timeStep / 2. * (vel_ + v0);

//    torque0_ = torque_;
//    torque_ = ibObj->torque();
//
//    omega_ += timeStep / (2 * ibObj->momentOfInertia()) * (torque0_ + torque_);

    ibObj->shape().move(pos_);
}
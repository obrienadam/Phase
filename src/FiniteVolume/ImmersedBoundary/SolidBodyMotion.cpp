#include "SolidBodyMotion.h"
#include "ImmersedBoundaryObject.h"

SolidBodyMotion::SolidBodyMotion(std::weak_ptr<ImmersedBoundaryObject> ibObj)
        :
        Motion(ibObj)
{
    if (ibObj.lock()->mass() <= 0.)
        throw Exception("SolidBodyMotion",
                        "SolidBodyMotion",
                        "must specify a non-zero density for immersed boundary object \"" + ibObj.lock()->name() +
                        "\".");

    force_ = Vector2D(0., 0.);
    torque_ = 0.;
}

void SolidBodyMotion::update(Scalar timeStep)
{
    auto ibObj = ibObj_.lock();

    //- Update translational motion
    Vector2D acc0 = acc_;
    force_ = ibObj->force();
    acc_ = force_ / ibObj->mass();

    Vector2D v0 = vel_;
    vel_ += timeStep / 2. * (acc_ + acc0);
    pos_ += timeStep / 2. * (vel_ + v0);

    //- Update rotational motion
    Scalar alpha0 = alpha_;
    torque_ = ibObj->torque();
    alpha_ = torque_ / ibObj->momentOfInertia();

    Scalar omega0 = omega_, theta0 = theta_;
    omega_ += timeStep / 2. * (alpha_ + alpha0);
    theta_ = std::fmod(theta_ + timeStep / 2. * (omega_ + omega0), 2. * M_PI);

    ibObj->shape().move(pos_);
    ibObj->shape().rotate(theta_ - theta0);
}
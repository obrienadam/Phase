#include "SolidBodyMotion.h"

SolidBodyMotion::SolidBodyMotion(
    std::weak_ptr<const ImmersedBoundaryObject> ibObj, const Vector2D &v0,
    bool constrainMotion, const Vector2D &motionAxis)
    : Motion(), ibObj_(ibObj), constrainMotion_(constrainMotion),
      motionAxis_(motionAxis.unitVec()) {
  if (ibObj.lock()->mass() <= 0.)
    throw Exception(
        "SolidBodyMotion", "SolidBodyMotion",
        "must specify a non-zero density for immersed boundary object \"" +
            ibObj.lock()->name() + "\".");

  pos_ = ibObj_.lock()->position();
  vel_ = constrainMotion_ ? v0 : dot(v0, motionAxis_) * motionAxis_;
  force_ = Vector2D(0., 0.);
  torque_ = 0.;
}

void SolidBodyMotion::update(Scalar timeStep) {
  auto ibObj = ibObj_.lock();

  //- Update translational motion
  Vector2D acc0 = acc_;

  force_ = constrainMotion_ ? dot(ibObj->force(), motionAxis_) * motionAxis_
                            : ibObj->force();

  acc_ = force_ / ibObj->mass();

  Vector2D v0 = vel_;
  vel_ += timeStep * (acc_ + acc0) / 2.;
  pos_ += timeStep * (vel_ + v0) / 2.;

  //- Update rotational motion
  //    Scalar alpha0 = alpha_;
  //    torque_ = ibObj->torque();
  //    alpha_ = torque_ / ibObj->momentOfInertia();

  //    Scalar omega0 = omega_, theta0 = theta_;
  //    omega_ += timeStep * (alpha_ + alpha0) / 2.;
  //    theta_ = std::fmod(theta_ + timeStep / 2. * (omega_ + omega0), 2. *
  //    M_PI);
}

void SolidBodyMotion::setMotionConstraint(const Vector2D &axis) {
  constrainMotion_ = true;
  motionAxis_ = axis.unitVec();
  acc_ = dot(acc_, motionAxis_) * motionAxis_;
  vel_ = dot(vel_, motionAxis_) * motionAxis_;
}

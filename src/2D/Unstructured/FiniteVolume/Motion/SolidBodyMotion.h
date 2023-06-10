#ifndef PHASE_SOLID_BODY_MOTION_H
#define PHASE_SOLID_BODY_MOTION_H

#include "FiniteVolume/ImmersedBoundary/ImmersedBoundaryObject.h"

#include "Motion.h"

class SolidBodyMotion : public Motion {
public:
  SolidBodyMotion(std::weak_ptr<const ImmersedBoundaryObject> ibObj,
                  const Vector2D &v0 = Vector2D(0., 0.),
                  bool constrainMotion = false,
                  const Vector2D &motionAxis = Vector2D(0., 1.));

  void update(Scalar timeStep);

  void setMotionConstraint(const Vector2D &axis);

private:
  Vector2D force_;

  Scalar torque_;

  bool constrainMotion_;

  Vector2D motionAxis_;

  std::weak_ptr<const ImmersedBoundaryObject> ibObj_;
};

#endif

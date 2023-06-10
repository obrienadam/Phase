#include "FiniteVolume/ImmersedBoundary/ImmersedBoundaryObject.h"

#include "OscillatingMotion.h"

OscillatingMotion::OscillatingMotion(const Point2D &pos, const Vector2D &freq,
                                     const Vector2D &amp, const Vector2D &phase,
                                     Scalar time)
    : Motion(pos), freq_(freq), amp_(amp), phase_(phase), time_(time) {
  pos0_ = pos;
}

void OscillatingMotion::update(Scalar timeStep) {
  time_ += timeStep;

  Vector2D omega(2 * M_PI * freq_.x, 2 * M_PI * freq_.y);

  pos_ = pos0_ + Vector2D(amp_.x * std::sin(omega.x * time_),
                          amp_.y * std::sin(omega.y * time_));
  vel_ = Vector2D(amp_.x * omega.x * std::cos(omega.x * time_),
                  amp_.y * omega.y * std::cos(omega.y * time_));
  acc_ = Vector2D(-amp_.x * std::pow(omega.x, 2) * std::sin(omega.x * time_),
                  -amp_.y * std::pow(omega.y, 2) * std::sin(omega.y * time_));
  // ibObj_.lock()->shape().move(x);
}

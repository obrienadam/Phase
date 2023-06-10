#include "MotionProfile.h"
#include "OscillatingMotion.h"
#include "TranslatingMotion.h"

MotionProfile::MotionProfile(const Point2D &pos, const Vector2D &vel,
                             const Vector2D &acc, Scalar theta, Scalar omega,
                             Scalar alpha)
    : Motion(pos, vel, acc, theta, omega, alpha), time_(0.) {}

void MotionProfile::update(Scalar timeStep) {
  auto it = std::upper_bound(
      timePoints_.begin(), timePoints_.end(), time_,
      [](Scalar time, const TimePoint &tp) { return time < tp.time; });

  if (it > timePoints_.begin() && !timePoints_.empty()) {
    --it;
    it->motion->init(pos_, vel_, acc_, theta_, omega_, alpha_);
    it->motion->update(timeStep);
    pos_ = it->motion->position();
    vel_ = it->motion->velocity();
    acc_ = it->motion->acceleration();
    theta_ = it->motion->theta();
    omega_ = it->motion->omega();
    alpha_ = it->motion->alpha();
  }
}

void MotionProfile::addMotion(Scalar startTime,
                              const std::shared_ptr<Motion> &motion) {
  timePoints_.push_back(TimePoint{startTime, motion});
  std::sort(timePoints_.begin(), timePoints_.end());
}

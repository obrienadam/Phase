#ifndef PHASE_MOTION_PROFILE_H
#define PHASE_MOTION_PROFILE_H

#include "Motion.h"

class MotionProfile: public Motion
{
public:

    MotionProfile(const Point2D &pos = Point2D(0., 0.),
                  const Vector2D &vel = Vector2D(0., 0.),
                  const Vector2D &acc = Vector2D(0., 0.),
                  Scalar theta = 0.,
                  Scalar omega = 0.,
                  Scalar alpha = 0.);

    virtual void update(Scalar timeStep) override;

    void addMotion(Scalar startTime, const std::shared_ptr<Motion> &motion);

protected:

    struct TimePoint
    {
        bool operator<(const TimePoint &rhs) const
        { return time < rhs.time; }

        Scalar time;

        std::shared_ptr<Motion> motion;
    };

    Scalar time_;

    std::vector<TimePoint> timePoints_;
};

#endif

#ifndef PHASE_OSCILLATING_MOTION_H
#define PHASE_OSCILLATING_MOTION_H

#include "Motion.h"

class OscillatingMotion: public Motion
{
public:

    OscillatingMotion(const Point2D& pos,
                      const Vector2D& freq,
                      const Vector2D& amp,
                      const Vector2D& phase = Vector2D(0, 0),
                      Scalar time = 0);

    void update(Scalar timeStep);

private:

    Point2D pos0_;

    Vector2D freq_, amp_, phase_;

    Scalar time_;
};

#endif

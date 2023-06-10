#ifndef TRANSLATING_MOTION_H
#define TRANSLATING_MOTION_H

#include "Motion.h"

class TranslatingMotion : public Motion {
public:
  TranslatingMotion(const Point2D &pos, const Vector2D &vel,
                    const Vector2D &acc = Vector2D(0., 0.));

  void update(Scalar timeStep);
};

#endif

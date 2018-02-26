#include "Motion.h"

Motion::Motion(const Point2D &pos,
               const Vector2D &vel,
               const Vector2D &acc,
               Scalar theta,
               Scalar omega,
               Scalar alpha)
        :
        pos_(pos),
        vel_(vel),
        acc_(acc),
        theta_(theta),
        omega_(omega),
        alpha_(alpha)
{

}
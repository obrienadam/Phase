#include "Motion.h"
#include "ImmersedBoundaryObject.h"

Motion::Motion(std::weak_ptr<ImmersedBoundaryObject> ibObj)
{
    ibObj_ = ibObj;
    acc_ = Vector2D(0., 0.);
    vel_ = Vector2D(0., 0.);
    pos_ = ibObj_.lock()->position();
}
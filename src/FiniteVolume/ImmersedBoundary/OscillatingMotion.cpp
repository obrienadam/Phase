#include "OscillatingMotion.h"

OscillatingMotion::OscillatingMotion(const Point2D &pos, const Vector2D &freq, const Vector2D &amp,
                                     const Vector2D &phase, Scalar time)
        :
        Motion(pos),
        freq_(freq),
        amp_(amp),
        phase_(phase),
        time_(time)
{

}

void OscillatingMotion::update(ImmersedBoundaryObject &ibObj, Scalar timeStep)
{
    time_ += timeStep;

    Vector2D omega(2 * M_PI * freq_.x, 2 * M_PI * freq_.y);

    Point2D x = pos_ + Vector2D(amp_.x * sin(omega.x * time_), amp_.y * sin(omega.y * time_));
    vel_ = Vector2D(amp_.x * omega.x * cos(omega.x * time_), amp_.y * omega.y * cos(omega.y * time_));
    acc_ = Vector2D(-amp_.x * pow(omega.x, 2) * sin(omega.x * time_), -amp_.y * pow(omega.y, 2) * sin(omega.y * time_));

    ibObj.shape().move(x);
}
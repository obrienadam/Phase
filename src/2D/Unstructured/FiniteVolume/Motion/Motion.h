#ifndef PHASE_MOTION_H
#define PHASE_MOTION_H

#include "Geometry/Point2D.h"

class Motion
{
public:

    Motion(const Point2D &pos = Point2D(0., 0.), const Vector2D &vel = Vector2D(0., 0.),
           const Vector2D &acc = Vector2D(0., 0.), Scalar theta = 0., Scalar omega = 0., Scalar alpha = 0.);

    virtual void update(Scalar timeStep) = 0;

    Point2D position() const
    { return pos_; }

    Vector2D acceleration(const Point2D &pt) const
    { return acc_ + alpha_ * (pt - pos_).tangentVec() + omega_ * omega_ * (pos_ - pt); }

    Vector2D velocity(const Point2D &pt) const
    { return vel_ + omega_ * (pt - pos_).tangentVec(); }

    Scalar alpha() const
    { return alpha_; }

    Scalar omega() const
    { return omega_; }

    Scalar theta() const
    { return theta_; }

protected:

    Scalar alpha_, omega_, theta_;
    Vector2D acc_, vel_, pos_;
};

#endif

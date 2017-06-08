#include "TranslatingMotion.h"

TranslatingMotion::TranslatingMotion(const Point2D &pos, const Vector2D &vel, const Vector2D &acc)
        :
        Motion(pos)
{
    vel_ = vel;
    acc_ = acc;
}

void TranslatingMotion::update(ImmersedBoundaryObject &ibObj, Scalar timeStep)
{
    pos_ += vel_*timeStep + acc_*timeStep*timeStep/2.;
    vel_ += acc_*timeStep;
    ibObj.shape().move(pos_);
}
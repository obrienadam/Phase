#include "TranslatingMotion.h"
#include "ImmersedBoundaryObject.h"

TranslatingMotion::TranslatingMotion(std::weak_ptr<ImmersedBoundaryObject> ibObj,
                                     const Vector2D &vel,
                                     const Vector2D &acc)
        :
        Motion(ibObj)
{
    vel_ = vel;
    acc_ = acc;
}

void TranslatingMotion::update(Scalar timeStep)
{
    pos_ += vel_*timeStep + acc_*timeStep*timeStep/2.;
    vel_ += acc_*timeStep;
    ibObj_.lock()->shape().move(pos_);
}
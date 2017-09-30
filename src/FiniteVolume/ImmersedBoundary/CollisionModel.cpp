#include "CollisionModel.h"

CollisionModel::CollisionModel(Scalar eps, Scalar range)
{
    eps_ = eps;
    range_ = range;
}

Vector2D CollisionModel::force(const ImmersedBoundaryObject &ibObjP, const ImmersedBoundaryObject &ibObjQ) const
{
    Vector2D xp = ibObjP.shape().nearestPoint(ibObjQ.shape());
    Vector2D xq = ibObjQ.shape().nearestPoint(ibObjP.shape());

    return (xp - xq) / eps_ * std::pow(std::max(0., -((xp - xq).mag() - range_)), 2);
}
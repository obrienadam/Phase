#include "BoundingBox.h"

BoundingBox::BoundingBox(const Point2D &lBound, const Point2D &uBound)
        :
        lBound_(lBound),
        uBound_(uBound)
{
    center_ = (lBound_ + uBound_) / 2.;
}

BoundingBox::Quadrant BoundingBox::getQuadrant(const Point2D &point) const
{
    Vector2D rVec = point - center_;

    if (rVec.x > 0. && rVec.y > 0.)
        return I;
    else if (rVec.x < 0. && rVec.y > 0.)
        return II;
    else if (rVec.x < 0. && rVec.y < 0.)
        return III;
    else
        return IV;
}

bool BoundingBox::isInBox(const Point2D &point) const
{
    return point.x < uBound_.x && point.y < uBound_.y && point.x > lBound_.x && point.y > lBound_.y;
}

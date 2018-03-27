#include "BoundingBox.h"

BoundingBox::BoundingBox(const Point2D &lBound, const Point2D &uBound)
        :
        lBound_(lBound),
        uBound_(uBound)
{
    center_ = (lBound_ + uBound_) / 2.;
}

BoundingBox::BoundingBox(const Point2D *points, size_t nPoints)
        :
        lBound_(points[0]),
        uBound_(points[0])
{
    for (size_t i = 0; i < nPoints; ++i)
    {
        const Point2D &point = points[i];

        if (point.x < lBound_.x)
            lBound_.x = point.x;
        if (point.y < lBound_.y)
            lBound_.y = point.y;

        if (point.x > uBound_.x)
            uBound_.x = point.x;
        if (point.y > uBound_.y)
            uBound_.y = point.y;
    }

    center_ = (lBound_ + uBound_) / 2.;
}

BoundingBox::Quadrant BoundingBox::getQuadrant(const Point2D &point) const
{
    Vector2D rVec = point - center_;

    if (rVec.x > 0. && rVec.y > 0.)
        return I;
    if (rVec.x < 0. && rVec.y > 0.)
        return II;
    if (rVec.x < 0. && rVec.y < 0.)
        return III;
    else
        return IV;
}

bool BoundingBox::isInBox(const Point2D &point) const
{
    return point.x < uBound_.x && point.y < uBound_.y && point.x > lBound_.x && point.y > lBound_.y;
}

std::string BoundingBox::toString() const
{
    using namespace std;

    return lBound_.toString() + " " + uBound_.toString();
}

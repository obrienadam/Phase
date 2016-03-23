#ifndef SHAPE_2D_H
#define SHAPE_2D_H

#include "Point2D.h"

class Shape2D
{
public:

    virtual Point2D centroid() const = 0;
    virtual Scalar area() const = 0;

    virtual bool isInside(const Point2D& testPoint) const = 0;
    bool isOutside(const Point2D& testPoint) const { return !isInside(testPoint); }

    virtual bool isOnEdge(const Point2D& testPoint) const = 0;
    virtual Point2D nearestIntersect(const Point2D& testPoint) const = 0;

    virtual void operator+=(const Vector2D& translationVec) = 0;
    virtual void operator-=(const Vector2D& translationVec) = 0;

private:
};

#endif

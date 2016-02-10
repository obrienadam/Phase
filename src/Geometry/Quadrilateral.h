#ifndef QUADRILATERAL_H
#define QUADRILATERAL_H

#include "Shape2D.h"

class Quadrilateral : public Shape2D
{
public:

    Quadrilateral() {}
    Quadrilateral(const Point2D& v1, const Point2D& v2, const Point2D& v3, const Point2D& v4);

    Scalar area() const;

    virtual bool isInside(const Point2D& testPoint) const;
    virtual bool isOnEdge(const Point2D& testPoint) const;

    virtual Point2D nearestIntersect(const Point2D& testPoint) const;

    virtual void operator+=(const Vector2D& translationVec);
    virtual void operator-=(const Vector2D& translationVec);

private:

    Point2D vertices_[4];
};

#endif

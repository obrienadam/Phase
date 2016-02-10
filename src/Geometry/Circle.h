#ifndef CIRCLE_H
#define CIRCLE_H

#include "Types.h"
#include "Shape2D.h"

class Circle : public Shape2D
{
public:

    Circle(Point2D center, Scalar radius);

    Scalar area() const;

    bool isInside(const Point2D &testPoint) const;
    bool isOnEdge(const Point2D &testPoint) const;
    Point2D nearestIntersect(const Point2D &testPoint) const;

    void operator+=(const Vector2D& translationVec);
    void operator-=(const Vector2D& translationVec);

private:

    Point2D center_;
    Scalar radius_;

};

#endif

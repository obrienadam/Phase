#ifndef CIRCLE_H
#define CIRCLE_H

#include<CGAL/Polygon_2.h>

#include "Types.h"
#include "Shape2D.h"
#include "Geometry.h"
#include "Polygon.h"

class Circle : public Shape2D
{
public:

    Circle(Point2D center, Scalar radius);

    void init(const Point2D &center, Scalar radius);

    const Point2D& centroid() const { return center_; }
    Scalar area() const;

    bool isInside(const Point2D &testPoint) const;
    bool isOnEdge(const Point2D &testPoint) const;
    Point2D nearestIntersect(const Point2D &testPoint) const;

    void operator+=(const Vector2D& translationVec);
    void operator-=(const Vector2D& translationVec);

    Kernel::Circle_2 cgalCircle() const { return Kernel::Circle_2(center_.cgalPoint(), radius_*radius_); }

private:

    Point2D center_;
    Scalar radius_, area_;

};

#endif

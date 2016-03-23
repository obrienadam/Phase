#ifndef POLYGON_H
#define POLYGON_H

#include <vector>

#include "Shape2D.h"

class Polygon : public Shape2D
{
public:

    Polygon(const std::vector<Point2D>& vertices);

    Point2D centroid() const { return centroid_; }
    Scalar area() const { return area_; }

    virtual bool isInside(const Point2D& testPoint) const;
    virtual bool isOnEdge(const Point2D& testPoint) const;

    virtual Point2D nearestIntersect(const Point2D& testPoint) const;

    virtual void operator+=(const Vector2D& translationVec);
    virtual void operator-=(const Vector2D& translationVec);

    bool isConvex() const;

private:

    Point2D centroid_;
    Scalar area_;
    std::vector<Point2D> vertices_;
};

#endif

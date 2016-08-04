#include <math.h>

#include "Circle.h"

Circle::Circle(Point2D center, Scalar radius)
    :
      center_(center),
      radius_(radius)
{
    area_ = M_PI*radius_*radius_;
}

void Circle::init(const Point2D &center, Scalar radius)
{
    center_ = center;
    radius_ = radius;
    area_ = M_PI*radius_*radius_;
}

Scalar Circle::area() const
{
    return area_;
}

void Circle::scale(Scalar factor)
{
    radius_ *= factor;
    area_ *= factor*factor;
}

bool Circle::isInside(const Point2D &testPoint) const
{
    return (testPoint - center_).magSqr() < radius_*radius_;
}

bool Circle::isOnEdge(const Point2D &testPoint) const
{
    return fabs((testPoint - center_).mag() - radius_) < 1e-12;
}

Point2D Circle::nearestIntersect(const Point2D &testPoint) const
{
    Vector2D rVec = testPoint - center_;
    return center_ + rVec.unitVec()*radius_;
}

std::pair<Point2D, bool> Circle::firstIntersect(Point2D ptA, Point2D ptB) const
{
    ptA -= center_;
    ptB -= center_;

    Scalar dx = ptB.x - ptA.x;
    Scalar dy = ptB.y - ptA.y;

    Scalar a = dx*dx + dy*dy;
    Scalar b = 2.*(ptA.x*dx + ptA.y*dy);
    Scalar c = ptA.x*ptA.x + ptA.y*ptA.y - radius_*radius_;
    Scalar delta = b*b - 4.*a*c;

    if(delta < 0.) //- no intersection
        return std::make_pair(Point2D(), false);

    Scalar t1 = (-b + sqrt(delta))/(2.*a),
            t2 = (-b - sqrt(delta))/(2.*a);

    ptA += center_;
    ptB += center_;

    return fabs(t1) < fabs(t2) ? std::make_pair(ptA + t1*(ptB - ptA), true) : std::make_pair(ptA + t2*(ptB - ptA), true);
}

boost::geometry::model::box<Point2D> Circle::boundingBox() const
{
    return boost::geometry::model::box<Point2D>(Point2D(center_.x - radius_,
                                                        center_.y - radius_),
                                                Point2D(center_.x + radius_,
                                                        center_.y + radius_));
}

void Circle::operator +=(const Vector2D& translationVec)
{
    center_ += translationVec;
}

void Circle::operator -=(const Vector2D& translationVec)
{
    center_ -= translationVec;
}

//- External functions

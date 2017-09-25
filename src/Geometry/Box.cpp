#include "Box.h"
#include "Polygon.h"

Box::Box(const Point2D &lower, const Point2D &upper)
        :
        lower_(lower),
        upper_(upper)
{
    if (upper_.x < lower_.x)
        std::swap(upper_.x, lower_.x);

    if (upper_.y < lower_.y)
        std::swap(upper_.y, lower_.y);

    Scalar w = upper_.x - lower_.x;
    Scalar h = upper_.y - lower_.y;

    centroid_ = (lower_ + upper_) / 2.;
    area_ = w * h;
    momentOfInertia_ = w * h / 12. * (w * w + h * h);
}

//- Tests
bool Box::isInside(const Point2D &point) const
{
    return point.x < upper_.x && point.x > lower_.x && point.y < upper_.y && point.y > lower_.y;
}

bool Box::isOnEdge(const Point2D &point) const
{
    return (point.x < upper_.x && point.x > lower_.x && (point.y == upper_.y || point.y == lower_.y)) ||
           (point.y < upper_.y && point.y > lower_.y && (point.x == upper_.x || point.x == lower_.x));
}

bool Box::isCovered(const Point2D &point) const
{
    return point.x <= upper_.x && point.x >= lower_.x && point.y <= upper_.y && point.y >= lower_.y;
}

//- Intersections
std::vector<Point2D> Box::intersections(const Line2D &line) const
{
    std::vector<Point2D> intersections;

    for (const LineSegment2D &edge: edges())
    {
        Vector2D r2 = edge.rVec();
        Vector2D b = edge.ptA() - line.r0();

        Scalar t = cross(line.d(), b) / cross(line.d(), -r2);

        if (t >= 0. && t <= 1.)
            intersections.push_back(edge.ptA() + t * r2);
    }

    return intersections;
}

std::vector<Point2D> Box::intersections(const LineSegment2D &line) const
{
    Vector2D r1 = line.rVec();
    std::vector<Point2D> intersections;

    for (const LineSegment2D &edge: edges())
    {
        Vector2D r2 = edge.rVec();
        Vector2D b = edge.ptA() - line.ptA();

        Scalar detA = cross(r1, -r2);

        Scalar t1 = cross(b, -r2) / detA;
        Scalar t2 = cross(r1, b) / detA;

        if ((t1 >= 0. && t1 <= 1.) && (t2 >= 0. && t2 <= 0.))
            intersections.push_back(line.ptA() + t1 * r1);
    }

    return intersections;
}

std::vector<Point2D> Box::intersections(const Ray2D& ray) const
{
    std::vector<Point2D> intersections;

    for (const LineSegment2D &edge: edges())
    {
        Scalar detA = cross(ray.r(), -edge.rVec());
        Vector2D b = edge.ptA() - ray.x0();

        Scalar t1 = cross(b, -edge.rVec()) / detA;
        Scalar t2 = cross(ray.r(), b) / detA;

        if (t1 > 0. && (t2 >= 0. && t2 <= 1.))
            intersections.push_back(ray(t1));
    }

    return intersections;
}

Point2D Box::nearestIntersect(const Point2D &point) const
{
    Point2D nearestXc;
    Scalar minDistSqr = std::numeric_limits<Scalar>::infinity();

    for (const LineSegment2D &edge: edges())
    {
        Vector2D rxc;

        if(edge.isBounded(point))
        {
            Vector2D r = point - edge.ptA();
            Vector2D t = edge.ptB() - edge.ptA();
            rxc = edge.ptA() + dot(r, t) * t / t.magSqr() - point;
        }
        else
        {
            rxc = ((edge.ptA() - point).magSqr() < (edge.ptB() - point).magSqr() ? edge.ptA(): edge.ptB()) - point;
        }

        Scalar distSqr = rxc.magSqr();

        if(distSqr < minDistSqr)
        {
            nearestXc = point + rxc;
            minDistSqr = distSqr;
        }
    }

    return nearestXc;
}

LineSegment2D Box::nearestEdge(const Point2D &point) const
{
    LineSegment2D nearestEdge;
    Scalar minDistSqr = std::numeric_limits<Scalar>::infinity();

    for (const LineSegment2D &edge: edges())
    {
        Scalar distSqr;

        if(edge.isBounded(point))
        {
            Vector2D r = point - edge.ptA();
            Vector2D t = edge.ptB() - edge.ptA();
            distSqr = (r - dot(r, t) * t / t.magSqr()).magSqr();
        }
        else
        {
            distSqr = std::min((edge.ptA() - point).magSqr(), (edge.ptB() - point).magSqr());
        }

        if(distSqr < minDistSqr)
        {
            nearestEdge = edge;
            minDistSqr = distSqr;
        }
    }

    return nearestEdge;
}

bool Box::intersects(const Shape2D &shape) const
{
    return isCovered(shape.nearestIntersect(centroid_));
}

void Box::scale(Scalar factor)
{
    upper_ += factor * (upper_ - centroid_);
    lower_ += factor * (lower_ - centroid_);

    Scalar w = upper_.x - lower_.x;
    Scalar h = upper_.y - lower_.y;

    area_ = w * h;
    momentOfInertia_ = w * h / 12. * (w * w + h * h);
}

void Box::rotate(Scalar theta)
{
    return;
}

Box &Box::move(const Point2D &pos)
{
    lower_ += (pos - centroid_);
    upper_ += (pos - centroid_);
    centroid_ = pos;
    return *this;
}

Box &Box::operator+=(const Vector2D &vec)
{
    centroid_ += vec;
    return *this;
}

Box &Box::operator-=(const Vector2D &vec)
{
    centroid_ -= vec;
    return *this;
}

boost::geometry::model::box<Point2D> Box::boundingBox() const
{
    return boost::geometry::model::box<Point2D>(
            lower_,
            upper_
    );
}

Polygon Box::polygonize() const
{
    return Polygon({
            lower_,
            Point2D(upper_.x, lower_.y),
            upper_,
            Point2D(lower_.x, upper_.y)
                                        });
}

std::array<Point2D, 4> Box::vertices() const
{
    return {
            lower_,
            Point2D(upper_.x, lower_.y),
            upper_,
            Point2D(lower_.x, upper_.y)
    };
}


std::array<LineSegment2D, 4> Box::edges() const
{
    return {
            LineSegment2D(lower_, Point2D(upper_.x, lower_.y)),
            LineSegment2D(Point2D(upper_.x, lower_.y), upper_),
            LineSegment2D(upper_, Point2D(lower_.x, upper_.y)),
            LineSegment2D(Point2D(lower_.x, upper_.y), lower_)
    };
}

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

    centroid_ = (lower_ + upper_) / 2.;
    area_ = (upper_.x - lower_.x) * (upper_.y - lower_.y);
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

bool Box::isBoundedBy(const Point2D &point, Scalar toler) const
{
    Point2D upper = Point2D(upper_.x + toler, upper_.y + toler);
    Point2D lower = Point2D(lower_.x - toler, lower_.y - toler);

    return point.x <= upper.x && point.x >= lower.x && point.y <= upper.y && point.y >= lower.y;
}

//- Intersections
std::vector<Point2D> Box::intersections(const Line2D &line) const
{
    std::vector<Point2D> intersections;

    for (const LineSegment2D &edge: edges())
    {
        Vector2D r2 = edge.ptB() - edge.ptA();
        Vector2D dr0 = edge.ptA() - line.r0();

        Scalar t2 = cross(line.d(), dr0) / cross(line.d(), -r2);

        if (edge.isBounded(edge.ptA() + t2 * r2))
            intersections.push_back(edge.ptA() + t2 * r2);
    }

    return intersections;
}

std::vector<Point2D> Box::intersections(const LineSegment2D &line) const
{
    Vector2D r1 = line.ptB() - line.ptA();
    std::vector<Point2D> intersections;

    for (const LineSegment2D &edge: edges())
    {
        Vector2D r2 = edge.ptB() - edge.ptA();
        Vector2D dr0 = edge.ptA() - line.ptA();

        Scalar t1 = cross(dr0, -r2) / cross(r1, -r2);
        Scalar t2 = cross(r1, dr0) / cross(r1, -r2);

        if (line.isBounded(line.ptA() + t1 * r1) && edge.isBounded(edge.ptA() + t2 * r2))
            intersections.push_back(line.ptA() + t1 * r1);
    }

    return intersections;
}

Point2D Box::nearestIntersect(const Point2D &point) const
{
    Point2D xc;
    Scalar minDistSqr = std::numeric_limits<Scalar>::infinity();

    for (const LineSegment2D &edge: edges())
    {
        if (!edge.isBounded(point))
            continue;

        Vector2D t = edge.ptB() - edge.ptA();
        t = dot(point - edge.ptA(), t) * t / t.magSqr();

        Vector2D n = point - edge.ptA() - t;

        if (n.magSqr() < minDistSqr)
        {
            minDistSqr = n.magSqr();
            xc = edge.ptA() + t;
        }
    }

    for (const Point2D &vertex: vertices())
        if ((point - vertex).magSqr() < minDistSqr)
        {
            minDistSqr = (point - vertex).magSqr();
            xc = vertex;
        }

    return xc;
}

LineSegment2D Box::nearestEdge(const Point2D &point) const
{
    LineSegment2D nearestEdge;
    Scalar minDistSqr = std::numeric_limits<Scalar>::infinity();

    for (const LineSegment2D &edge: edges())
    {
        if (!edge.isBounded(point))
            continue;

        Vector2D t = edge.ptB() - edge.ptA();
        t = dot(point - edge.ptA(), t) * t / t.magSqr();

        Vector2D n = point - edge.ptA() - t;

        if (n.magSqr() < minDistSqr)
        {
            nearestEdge = edge;
            minDistSqr = n.magSqr();
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

    area_ = (upper_.x - lower_.x) * (upper_.y - lower_.y);
}

void Box::rotate(Scalar theta)
{
    return;
}

Shape2D &Box::operator+=(const Vector2D &vec)
{
    centroid_ += vec;
    return *this;
}

Shape2D &Box::operator-=(const Vector2D &vec)
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
    return Polygon(std::vector<Point2D>({
            lower_,
            Point2D(upper_.x, lower_.y),
            upper_,
            Point2D(lower_.x, upper_.y)
                                        }));
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

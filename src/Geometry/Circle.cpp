#include "Circle.h"

Circle::Circle(Point2D center, Scalar radius)
        :
        center_(center),
        radius_(radius)
{
    area_ = M_PI * radius_ * radius_;
}

void Circle::init(const Point2D &center, Scalar radius)
{
    center_ = center;
    radius_ = radius;
    area_ = M_PI * radius_ * radius_;
}

//- Tests
bool Circle::isInside(const Point2D &testPoint) const
{
    return (testPoint - center_).magSqr() < radius_ * radius_;
}

bool Circle::isOnEdge(const Point2D &point) const
{
    return (point - center_).magSqr() == radius_ * radius_;
}

bool Circle::isCovered(const Point2D &point) const
{
    return (point - center_).magSqr() <= radius_ * radius_;
}

//- Intersections
std::vector<Point2D> Circle::intersections(const Line2D &line) const
{
    Vector2D r = line.r0() - center_;

    Scalar a = line.d().magSqr();
    Scalar b = 2*dot(r, line.d());
    Scalar c = r.magSqr() - radius_*radius_;
    Scalar disc = b*b - 4*a*c;

    if(disc < 0.)
        return std::vector<Point2D>();
    else if(disc == 0.)
        return std::vector<Point2D>({line(-b/(2.*a))});

    Scalar t1 = (-b - sqrt(disc))/(2*a);
    Scalar t2 = (-b + sqrt(disc))/(2*a);

    return std::vector<Point2D>({line(t1), line(t2)});
}

std::vector<Point2D> Circle::intersections(const LineSegment2D &line) const
{
    std::vector<Point2D> xc;

    for(const Point2D& pt: intersections(Line2D(line.ptA(), (line.ptB() - line.ptA()).normalVec())))
        if(line.isBounded(pt))
            xc.push_back(pt);

    return xc;
}

Point2D Circle::nearestIntersect(const Point2D &point) const
{
    return (point - center_).unitVec() * radius_ + center_;
}

LineSegment2D Circle::nearestEdge(const Point2D &point) const
{
    const Point2D xc = nearestIntersect(point);
    const Vector2D utan = (center_ - xc).tangentVec().unitVec();

    return LineSegment2D(
            xc + 0.5 * utan,
            xc - 0.5 * utan
    );
}

bool Circle::intersects(const Shape2D &shape) const
{
    return isCovered(shape.nearestIntersect(center_));
}

//- Transformations
void Circle::scale(Scalar factor)
{
    radius_ *= factor;
    area_ *= factor * factor;
}

//- Translations
Circle &Circle::move(const Point2D& pos)
{
    center_ = pos;
    return *this;
}

Circle &Circle::operator+=(const Vector2D &translationVec)
{
    center_ += translationVec;
    return *this;
}

Circle &Circle::operator-=(const Vector2D &translationVec)
{
    center_ -= translationVec;
    return *this;
}

boost::geometry::model::box<Point2D> Circle::boundingBox() const
{
    return boost::geometry::model::box<Point2D>(Point2D(center_.x - radius_,
                                                        center_.y - radius_),
                                                Point2D(center_.x + radius_,
                                                        center_.y + radius_));
}

Polygon Circle::polygonize(int nVerts) const
{
    std::vector<Point2D> verts(nVerts);

    for (int i = 0; i < nVerts; ++i)
    {
        const Scalar theta = 2. * M_PI * i / Scalar(nVerts);
        verts[i] = center_ + Point2D(radius_ * cos(theta), radius_ * sin(theta));
    }

    return Polygon(verts);
}

Polygon Circle::polygonize() const
{
    return polygonize(200);
}

//- External functions

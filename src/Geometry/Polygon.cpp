#include <math.h>

#include "Polygon.h"
#include "Exception.h"

Polygon::Polygon()
{
    area_ = 0.;
    centroid_ = Point2D(0., 0.);
}

Polygon::Polygon(const std::initializer_list<Point2D> &vertices)
    :
      Polygon(vertices.begin(), vertices.end())
{

}

Polygon::Polygon(const boost::geometry::model::ring<Point2D, false, true> &boostRing)
        :
        poly_(boostRing)
{
    init();
}

Scalar Polygon::perimeter() const
{
    return boost::geometry::perimeter(poly_);
}

//- Tests
bool Polygon::isInside(const Point2D &testPoint) const
{
    return boost::geometry::within(testPoint, poly_);
}

bool Polygon::isOnEdge(const Point2D &testPoint) const
{
    return boost::geometry::covered_by(testPoint, poly_) && !isInside(testPoint);
}

bool Polygon::isCovered(const Point2D &point) const
{
    return boost::geometry::covered_by(point, poly_);
}

bool Polygon::isValid() const
{
    return boost::geometry::is_valid(poly_);
}

bool Polygon::isEmpty() const
{
    return vertices().size() == 0;
}

//- Intersections
std::vector<Point2D> Polygon::intersections(const Line2D &line) const
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

std::vector<Point2D> Polygon::intersections(const LineSegment2D &line) const
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

Point2D Polygon::nearestIntersect(const Point2D &point) const
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

LineSegment2D Polygon::nearestEdge(const Point2D &point) const
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

bool Polygon::intersects(const Shape2D &shape) const
{
    return boost::geometry::intersects(poly_, shape.polygonize().boostRing());
}

//- Transformations
void Polygon::scale(Scalar factor)
{
    for (Point2D &vtx: poly_)
        vtx = factor * (vtx - centroid_) + centroid_;

    init();
}

void Polygon::rotate(Scalar theta)
{
    for (Point2D &vtx: poly_)
        vtx = (vtx - centroid_).rotate(theta) + centroid_;
}

Polygon Polygon::scale(Scalar factor) const
{
    std::vector<Point2D> verts;
    verts.reserve(vertices().size());

    for (const Point2D &vtx: poly_)
        verts.push_back(factor * (vtx - centroid_) + centroid_);

    return Polygon(verts.begin(), verts.end());
}

//- Translations
Polygon &Polygon::move(const Point2D& pos)
{
    for (Point2D &vtx: poly_)
        vtx += (pos - centroid_);

    centroid_ = pos;

    return *this;
}

Polygon &Polygon::operator+=(const Vector2D &translationVec)
{
    for (Point2D &vtx: poly_)
        vtx += translationVec;

    centroid_ += translationVec;

    return *this;
}

Polygon &Polygon::operator-=(const Vector2D &translationVec)
{
    return operator+=(-translationVec);
}

//- Bounding box
boost::geometry::model::box<Point2D> Polygon::boundingBox() const
{
    boost::geometry::model::box<Point2D> box;
    boost::geometry::envelope(poly_, box);

    return box;
}

std::vector<LineSegment2D> Polygon::edges() const
{
    std::vector<LineSegment2D> edges;

    auto vtxA = poly_.begin();
    auto vtxB = vtxA + 1;

    for (; vtxB != vertices().end(); ++vtxA, ++vtxB)
        edges.push_back(LineSegment2D(*vtxA, *vtxB));

    return edges;
}

//- Protected methods

void Polygon::init()
{
    if (poly_.size() > 0)
    {
        boost::geometry::unique(poly_);
        boost::geometry::correct(poly_);
        area_ = boost::geometry::area(poly_);
        boost::geometry::centroid(poly_, centroid_);
    }
    else
    {
        area_ = 0.;
        centroid_ = Point2D(0., 0.);
    }

    valid_ = isValid();
    simple_ = boost::geometry::is_simple(poly_);
}

//- External functions
Polygon intersectionPolygon(const Polygon &pgnA, const Polygon &pgnB)
{
    std::vector<boost::geometry::model::ring<Point2D, false, true> > pgn;

    boost::geometry::intersection(pgnA.boostRing(), pgnB.boostRing(), pgn);

    if (pgn.size() == 0)
        return Polygon();
    else if (pgn.size() > 1)
        throw Exception("Polygon", "intersectionPolygon", "there are two polygons!");

    return Polygon(pgn.front());
}

Polygon difference(const Polygon &pgnA, const Polygon &pgnB)
{
    std::vector<boost::geometry::model::ring<Point2D, false, true> > pgn;
    boost::geometry::difference(pgnA.boostRing(), pgnB.boostRing(), pgn);

    if (pgn.size() == 0)
        return Polygon();
    else if (pgn.size() > 1)
        throw Exception("Polygon", "difference", "there is more than one output polygon!");

    return Polygon(pgn.front());
}

#include <math.h>

#include "Polygon.h"

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
        Vector2D r2 = edge.rVec();
        Vector2D b = edge.ptA() - line.r0();

        Scalar t = cross(line.d(), b) / cross(line.d(), -r2);

        if (t >= 0. && t <= 1.)
            intersections.push_back(edge.ptA() + t * r2);
    }

    return intersections;
}

std::vector<Point2D> Polygon::intersections(const LineSegment2D &line) const
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

std::vector<Point2D> Polygon::intersections(const Ray2D &ray) const
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

Point2D Polygon::nearestIntersect(const Point2D &point) const
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

LineSegment2D Polygon::nearestEdge(const Point2D &point) const
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
std::vector<Polygon> intersection(const Polygon &pgnA, const Polygon &pgnB)
{
    std::vector<boost::geometry::model::ring<Point2D, false, true> > boostPolygons;
    boost::geometry::intersection(pgnA.boostRing(), pgnB.boostRing(), boostPolygons);
    std::vector<Polygon> polygons(boostPolygons.size());

    std::transform(boostPolygons.begin(), boostPolygons.end(), polygons.begin(),
                   [](const boost::geometry::model::ring<Point2D, false, true>& ring) {
                       return Polygon(ring);
                   });

    return polygons;
}

std::vector<Polygon> difference(const Polygon &pgnA, const Polygon &pgnB)
{
    std::vector<boost::geometry::model::ring<Point2D, false, true> > boostPolygons;
    boost::geometry::difference(pgnA.boostRing(), pgnB.boostRing(), boostPolygons);
    std::vector<Polygon> polygons(boostPolygons.size());

    std::transform(boostPolygons.begin(), boostPolygons.end(), polygons.begin(),
                   [](const boost::geometry::model::ring<Point2D, false, true>& ring) {
                       return Polygon(ring);
                   });

    return polygons;
}

#include <math.h>

#include "Polygon.h"
#include "Exception.h"

Polygon::Polygon()
{
    area_ = 0.;
    centroid_ = Point2D(0., 0.);
}

Polygon::Polygon(const std::vector<Point2D> &vertices)
{
    for(const Point2D& vtx: vertices)
        boost::geometry::append(poly_, vtx);

    init();
}

Polygon::Polygon(const boost::geometry::model::polygon<Point2D, false, true> &boostPgn)
    :
      poly_(boostPgn)
{
    init();
}

//- Tests
bool Polygon::isInside(const Point2D& testPoint) const
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
std::vector<Point2D> Polygon::intersections(const Line2D& line) const
{
    auto vtx0 = vertices().begin();
    auto vtx1 = vtx0 + 1;
    std::vector<Point2D> intersections;

    for(; vtx1 != vertices().end(); ++vtx1, ++vtx0)
    {
        //- check for vertex intersections
        if(line.isApproximatelyOnLine(*vtx1))
            intersections.push_back(*vtx1);
        else if(line.isApproximatelyOnLine(*vtx0))
            continue;
        else // check for edge intersection
        {
            const Line2D edge(*vtx0, (*vtx1 - *vtx0).tangentVec());
            const Point2D xc = Line2D::intersection(edge, line).first;

            if((xc - *vtx0).magSqr() < (*vtx1 - *vtx0).magSqr() && dot(*vtx1 - *vtx0, xc - *vtx0) > 0.)
                intersections.push_back(xc);
        }
    }

    return intersections;
}

Point2D Polygon::nearestIntersect(const Point2D& point) const
{
    auto vtx0 = vertices().begin();
    auto vtx1 = vtx0 + 1;
    Scalar minDistSqr = (*vtx0 - point).magSqr();
    Vector2D nearestPoint = *vtx0;

    for(; vtx1 != vertices().end(); ++vtx0, ++vtx1)
    {
        auto edge = std::make_pair(*vtx0, *vtx1);
        const Vector2D edgeVec = edge.second - edge.first;

        Vector2D rel = point - edge.first;
        Vector2D tan = dot(rel, edgeVec)*edgeVec/edgeVec.magSqr();
        Vector2D norm = rel - tan;

        Scalar distSqr = norm.magSqr();
        const Point2D xc = edge.first + tan;

        const Scalar t = dot(xc - edge.first, tan.unitVec());
        const Scalar tmax = edgeVec.mag();

        if(distSqr < minDistSqr && t >= 0 && t < tmax)
        {
            minDistSqr = distSqr;
            nearestPoint = edge.first + tan;
        }

        distSqr = (edge.first - point).magSqr();

        if(distSqr < minDistSqr)
        {
            minDistSqr = distSqr;
            nearestPoint = edge.first;
        }
    }

    return nearestPoint;
}

std::pair<Point2D, Point2D> Polygon::nearestEdge(const Point2D& point) const
{
    auto vtx0 = vertices().begin();
    auto vtx1 = vtx0 + 1;
    Scalar minDistSqr = (*vtx0 - point).magSqr();
    std::pair<Point2D, Point2D> nearestEdge = std::make_pair(*vtx0, *vtx1);

    for(; vtx1 != vertices().end(); ++vtx0, ++vtx1)
    {
        auto edge = std::make_pair(*vtx0, *vtx1);
        const Vector2D edgeVec = edge.second - edge.first;

        Vector2D rel = point - edge.first;
        Vector2D tan = dot(rel, edgeVec)*edgeVec/edgeVec.magSqr();
        Vector2D norm = rel - tan;

        Scalar distSqr = norm.magSqr();
        const Point2D xc = edge.first + tan;

        const Scalar t = dot(xc - edge.first, tan.unitVec());
        const Scalar tmax = edgeVec.mag();

        if(distSqr < minDistSqr && t >= 0 && t < tmax)
        {
            minDistSqr = distSqr;
            nearestEdge = edge;
        }

        distSqr = (edge.first - point).magSqr();

        if(distSqr < minDistSqr)
        {
            minDistSqr = distSqr;
            nearestEdge = edge;
        }
    }

    return nearestEdge;
}

//- Transformations
void Polygon::scale(Scalar factor)
{
    for(Point2D &vtx: boost::geometry::exterior_ring(poly_))
        vtx = factor*(vtx - centroid_) + centroid_;

    init();
}

void Polygon::rotate(Scalar theta)
{
    for(Point2D &vtx: boost::geometry::exterior_ring(poly_))
        vtx = (vtx - centroid_).rotate(theta) + centroid_;
}

Polygon Polygon::scale(Scalar factor) const
{
    std::vector<Point2D> verts;
    verts.reserve(vertices().size());

    for(const Point2D &vtx: boost::geometry::exterior_ring(poly_))
        verts.push_back(factor*(vtx - centroid_) + centroid_);

    return Polygon(verts);
}

//- Translations
void Polygon::operator+=(const Vector2D& translationVec)
{
    for(Point2D &vtx: boost::geometry::exterior_ring(poly_))
        vtx += translationVec;

    centroid_ += translationVec;
}

void Polygon::operator-=(const Vector2D& translationVec)
{
    operator +=(-translationVec);
}

//- Bounding box
boost::geometry::model::box<Point2D> Polygon::boundingBox() const
{
    boost::geometry::model::box<Point2D> box;
    boost::geometry::envelope(boostPolygon(), box);

    return box;
}

//- Protected methods

void Polygon::init()
{
    if(boost::geometry::exterior_ring(poly_).size() > 0)
    {
        boost::geometry::correct(poly_);
        area_ = boost::geometry::area(poly_);
        boost::geometry::centroid(poly_, centroid_);
    }
    else
    {
        area_ = 0.;
        centroid_ = Point2D(0., 0.);
    }
}

//- External functions
Polygon intersectionPolygon(const Polygon &pgnA, const Polygon &pgnB)
{
    std::vector< boost::geometry::model::polygon<Point2D, false, true> > pgn;

    boost::geometry::intersection(pgnA.boostPolygon(), pgnB.boostPolygon(), pgn);

    if(pgn.size() == 0)
        return Polygon();
    else if (pgn.size() > 1)
        throw Exception("Polygon", "intersectionPolygon", "there are two polygons!");

    return Polygon(pgn.front());
}

Polygon clipPolygon(const Polygon& pgn, const Line2D& line)
{
    std::vector<Point2D> verts;

    for(auto vertIt = pgn.vertices().begin(); vertIt != pgn.vertices().end() - 1; ++vertIt)
    {
        const Vector2D& vtx = *vertIt;
        const Vector2D& nextVtx = *(vertIt + 1);
        Line2D edgeLine = Line2D(vtx, (nextVtx - vtx).normalVec());

        if(line.isApproximatelyOnLine(vtx))
        {
            verts.push_back(vtx); // special case
            continue;
        }
        else if(line.isBelowLine(vtx))
            verts.push_back(vtx);

        std::pair<Point2D, bool> xc = Line2D::intersection(line, edgeLine);

        if(xc.second) // the lines are not paralell, ie xc is valid
        {
            Scalar l = (nextVtx - vtx).magSqr();
            Scalar x = dot(nextVtx - vtx, xc.first - vtx);

            if(x < l && x > 0 && !(xc.first == vtx || xc.first == nextVtx)) // intersection is on the segment
                verts.push_back(xc.first);
        }
    }

    return Polygon(verts);
}

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

bool Polygon::isInside(const Point2D& testPoint) const
{
    return boost::geometry::within(testPoint, poly_);
}

bool Polygon::isOnEdge(const Point2D &testPoint) const
{
    return boost::geometry::covered_by(testPoint, poly_) && !isInside(testPoint);
}

bool Polygon::isValid() const
{
    return boost::geometry::is_valid(poly_);
}

bool Polygon::isEmpty() const
{
    return vertices().size() == 0;
}

Point2D Polygon::nearestIntersect(const Point2D& testPoint) const
{
    auto vtxIter = begin();
    Point2D vtx0 = *(vtxIter++);

    Point2D currNearestIntersect(INFINITY, INFINITY);
    for(; vtxIter != end(); vtxIter++)
    {
        Point2D vtx1 = *(vtxIter);

        Point2D ptVec = testPoint - vtx0;
        Point2D segmentVec = vtx1 - vtx0;

        Scalar l = dot(ptVec, segmentVec.unitVec());

        Point2D intersect;

        if(l <= 0)
            intersect = vtx0;
        else if (l*l > segmentVec.magSqr())
            intersect = vtx1;
        else
            intersect = vtx0 + l*segmentVec.unitVec();

        currNearestIntersect = (intersect - testPoint).magSqr() < (currNearestIntersect - testPoint).magSqr() ? intersect : currNearestIntersect;

        vtx0 = vtx1;
    }

    return currNearestIntersect;
}

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

Polygon Polygon::scale(Scalar factor) const
{
    std::vector<Point2D> verts;
    verts.reserve(vertices().size());

    for(const Point2D &vtx: boost::geometry::exterior_ring(poly_))
        verts.push_back(factor*(vtx - centroid_) + centroid_);

    return Polygon(verts);
}

void Polygon::rotate(Scalar theta)
{
    for(Point2D &vtx: boost::geometry::exterior_ring(poly_))
        vtx = (vtx - centroid_).rotate(theta) + centroid_;
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

    for(auto vertIt = pgn.begin(); vertIt != pgn.end() - 1; ++vertIt)
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

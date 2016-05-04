#include <math.h>

#include "Polygon.h"
#include "Exception.h"

Polygon::Polygon()
{

}

Polygon::Polygon(const std::vector<Point2D> &vertices)
{
    for(const Point2D& vtx: vertices)
        boost::geometry::append(poly_, vtx);

    boost::geometry::append(poly_, vertices.front()); // weird that this is necessary

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

void Polygon::scale(Scalar factor)
{
    for(Point2D &vtx: boost::geometry::exterior_ring(poly_))
        vtx = factor*(vtx - centroid_) + centroid_;

    area_ = boost::geometry::area(poly_);
}

void Polygon::rotate(Scalar theta)
{
    for(Point2D &vtx: boost::geometry::exterior_ring(poly_))
        vtx = (vtx - centroid_).rotate(theta) + centroid_;
}

//- Protected methods

void Polygon::init()
{
    area_ = boost::geometry::area(poly_);
    boost::geometry::centroid(poly_, centroid_);
}

//- External functions

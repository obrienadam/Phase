#include <math.h>

#include "Polygon.h"

Polygon::Polygon(const std::vector<Point2D> &vertices)
{
    vertices_ = vertices;
    const int nPts = vertices_.size();

    area_ = 0.;

    for(int i = 0; i < nPts; ++i)
        area_ += vertices_[i].x*vertices_[(i + 1)%nPts].y - vertices_[(i + 1)%nPts].x*vertices_[i].y;

    area_ *= 0.5;

    Scalar coeff = 1./(6.*area_), ai;
    centroid_ = Point2D(0., 0.);

    for(int i = 0; i < nPts; ++i)
    {
        ai = vertices_[i].x*vertices_[(i + 1)%nPts].y - vertices_[(i + 1)%nPts].x*vertices_[i].y;
        centroid_.x += (vertices_[i].x + vertices_[(i + 1)%nPts].x)*ai;
        centroid_.y += (vertices_[i].y + vertices_[(i + 1)%nPts].y)*ai;
    }

    centroid_ *= coeff;
}

bool Polygon::isInside(const Point2D& testPoint) const
{

}

bool Polygon::isOnEdge(const Point2D& testPoint) const
{

}

Point2D Polygon::nearestIntersect(const Point2D& testPoint) const
{

}

void Polygon::operator+=(const Vector2D& translationVec)
{
    centroid_ += translationVec;

    for(Point2D& pt: vertices_)
        pt += translationVec;
}

void Polygon::operator-=(const Vector2D& translationVec)
{
    operator +=(-translationVec);
}

bool Polygon::isConvex() const
{
    Vector2D prevVec = vertices_[1] - vertices_[0];

    for(int i = 1, nPts = vertices_.size(); i < nPts; ++i)
    {
        Vector2D currVec = vertices_[(i + 1)%nPts] - vertices_[i];

        if(cross(prevVec, currVec) < 0.)
            return false;

        prevVec = currVec;
    }

    return true;
}

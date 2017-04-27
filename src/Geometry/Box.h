#ifndef BOX_H
#define BOX_H

#include "Shape2D.h"

class Box : public Shape2D
{
public:
    Box(const Point2D &lower = Point2D(), const Point2D &upper = Point2D());

    Type type() const
    { return BOX; }

    //- Shape parameters
    const Point2D &centroid() const
    { return centroid_; }

    Scalar area() const
    { return area_; }

    //- Tests
    bool isInside(const Point2D &point) const;

    bool isOnEdge(const Point2D &point) const;

    bool isCovered(const Point2D &point) const;

    bool isBoundedBy(const Point2D &point, Scalar toler) const;

    //- Intersections
    std::vector<Point2D> intersections(const Line2D &line) const;

    std::vector<Point2D> intersections(const LineSegment2D &line) const;

    Point2D nearestIntersect(const Point2D &point) const;

    LineSegment2D nearestEdge(const Point2D &point) const;

    bool intersects(const Shape2D &shape) const;

    //- Transformations
    void scale(Scalar factor);

    void rotate(Scalar theta);

    //- Translations
    Shape2D &operator+=(const Vector2D &vec);

    Shape2D &operator-=(const Vector2D &vec);

    //- Bounding box
    boost::geometry::model::box<Point2D> boundingBox() const;

    //- Convenience
    Polygon polygonize() const;

    //- Misc
    std::array<Point2D, 4> vertices() const;

    std::array<LineSegment2D, 4> edges() const;

    const Vector2D &lower() const
    { return lower_; }

    const Vector2D &upper() const
    { return upper_; }

private:
    Point2D centroid_, lower_, upper_;
    Scalar area_;
};


#endif

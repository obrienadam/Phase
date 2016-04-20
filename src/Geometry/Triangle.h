#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Geometry.h"
#include "Shape2D.h"

class Triangle : public Shape2D
{
public:

    Triangle(const Point2D vertices[]);
    Triangle(const Point2D& vtx1, const Point2D& vtx2, const Point2D& vtx3);
    Triangle(const Kernel::Triangle_2& cgalTriangle);

    Kernel::Triangle_2 cgalTriangle() const;

    Scalar area() const { return area_; }
    const Point2D& centroid() const { return centroid_; }

private:

    void init();

    Point2D vertices_[3];
    Scalar area_;
    Point2D centroid_;
};

#endif

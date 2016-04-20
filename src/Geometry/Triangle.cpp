#include "Triangle.h"

Triangle::Triangle(const Point2D vertices[])
{
    for(int i = 0; i < 3; ++i)
        vertices_[i] = vertices[i];

    init();
}

Triangle::Triangle(const Point2D &vtx1, const Point2D &vtx2, const Point2D &vtx3)
{
    vertices_[0] = vtx1;
    vertices_[1] = vtx2;
    vertices_[2] = vtx3;

    init();
}

Triangle::Triangle(const Kernel::Triangle_2 &cgalTriangle)
    :
      Triangle(
          cgalTriangle.vertex(0),
          cgalTriangle.vertex(1),
          cgalTriangle.vertex(2)
          )
{

}

Kernel::Triangle_2 Triangle::cgalTriangle() const
{
    return Kernel::Triangle_2(
                vertices_[0].cgalPoint(),
            vertices_[1].cgalPoint(),
            vertices_[2].cgalPoint()
            );
}

//- Private methods

void Triangle::init()
{
    area_ = cross(
                vertices_[1] - vertices_[0],
            vertices_[2] - vertices_[0]
                )/2.;

    centroid_ = (vertices_[0] + vertices_[1] + vertices_[2])/3.;
}

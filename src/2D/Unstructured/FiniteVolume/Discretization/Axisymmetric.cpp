#include "Axisymmetric.h"

Point2D axi::map(const Point2D &pt, const Vector2D &zaxis)
{
    return Point2D(-cross(pt, zaxis), dot(pt, zaxis));
}

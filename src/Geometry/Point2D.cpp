#include "Point2D.h"
#include "Vector2D.h"

Point2D::Point2D(const Vector2D &other)
        :
        Kernel::Point_2(other.x(), other.y())
{

}

std::string Point2D::toString() const
{
    return "(" + std::to_string(x()) + "," + std::to_string(y()) + ")";
}

//- External functions

Point2D operator+(const Point2D &lhs, const Point2D &rhs)
{
    return Point2D(lhs.x() + rhs.x(), lhs.y() + rhs.y());
}

Point2D operator-(const Point2D &lhs, const Point2D &rhs)
{
    return Point2D(lhs.x() - rhs.x(), lhs.y() - rhs.y());
}

Point2D operator*(const Scalar &lhs, const Point2D &rhs)
{
    return Point2D(lhs * rhs.x(), lhs * rhs.y());
}

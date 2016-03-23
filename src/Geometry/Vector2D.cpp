#include <math.h>

#include "Vector2D.h"
#include "Exception.h"

Vector2D::Vector2D(Scalar x, Scalar y)
    :
      x(x),
      y(y)
{

}

Scalar Vector2D::mag() const
{
    return sqrt(x*x + y*y);
}

Scalar Vector2D::magSqr() const
{
    return x*x + y*y;
}

Vector2D Vector2D::unitVec() const
{
    Scalar invMag = 1./mag();
    return Vector2D(x*invMag, y*invMag);
}

Vector2D Vector2D::normalVec() const
{
    return Vector2D(y, -x);
}

Scalar Vector2D::angle() const
{
    return atan2(y, x);
}

Scalar Vector2D::angle(const Vector2D& other) const
{
    return atan2(y - other.y, x - other.x);
}

Vector2D Vector2D::rotate(Scalar theta) const
{
    Scalar cosTheta = cos(theta), sinTheta = sin(theta);
    return Vector2D(x*cosTheta - y*sinTheta, x*sinTheta + y*cosTheta);
}

Vector2D Vector2D::transform(const Vector2D &uPrime) const
{
    return rotate((*this - uPrime).angle());
}

std::string Vector2D::toString() const
{
    using namespace std;

    return "(" + to_string(x) + ", " + to_string(y) + ")";
}

//- Operators

Scalar& Vector2D::operator ()(int component)
{
    switch(component)
    {
    case 0: return x;
    case 1: return y;
    default: throw Exception("Vector2D", "operator()", "invalid component.");
    }
}

Scalar Vector2D::operator ()(int component) const
{
    switch(component)
    {
    case 0: return x;
    case 1: return y;
    default: throw Exception("Vector2D", "operator()", "invalid component.");
    }
}

Vector2D& Vector2D::operator +=(const Vector2D& other)
{
    x += other.x;
    y += other.y;
    return *this;
}

Vector2D& Vector2D::operator -=(const Vector2D& other)
{
    x -= other.x;
    y -= other.y;
    return *this;
}

Vector2D& Vector2D::operator *=(Scalar other)
{
    x *= other;
    y *= other;
    return *this;
}

Vector2D& Vector2D::operator /=(Scalar other)
{
    x /= other;
    y /= other;
    return *this;
}

bool Vector2D::operator <(const Vector2D& rhs)
{
    return x < rhs.x || (x == rhs.x && y < rhs.y) ? true : false;
}

//- External functions

std::ostream& operator<<(std::ostream& os, const Vector2D& vec)
{
    return os << "(" << vec.x << ", " << vec.y << ")";
}

Vector2D operator+(Vector2D lhs, const Vector2D& rhs)
{
    lhs += rhs;
    return lhs;
}

Vector2D operator-(Vector2D lhs, const Vector2D& rhs)
{
    lhs -= rhs;
    return lhs;
}

Vector2D operator-(const Vector2D& rhs)
{
    return Vector2D(-rhs.x, -rhs.y);
}

Vector2D operator*(Vector2D lhs, Scalar rhs)
{
    lhs *= rhs;
    return lhs;
}

Vector2D operator*(Scalar lhs, Vector2D rhs)
{
    return rhs*lhs;
}

Vector2D operator/(Vector2D lhs, Scalar rhs)
{
    lhs /= rhs;
    return lhs;
}

Scalar dot(const Vector2D &u, const Vector2D &v)
{
    return u.x*v.x + u.y*v.y;
}

Scalar cross(const Vector2D &u, const Vector2D &v)
{
    return u.x*v.y - u.y*v.x;
}

std::string to_string(const Vector2D &vec)
{
    return "(" + std::to_string(vec.x) + ", " + std::to_string(vec.y) + ")";
}



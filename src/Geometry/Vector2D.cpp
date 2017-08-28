#include <math.h>

#include <boost/algorithm/string.hpp>

#include "Vector2D.h"
#include "Exception.h"

const Scalar Vector2D::EPSILON_ =
        10. * std::numeric_limits<Scalar>::epsilon(); // This parameter governs the precision of vector operations

Vector2D::Vector2D(std::string vecStr)
{
    using namespace std;
    using namespace boost::algorithm;

    vecStr = vecStr.substr(vecStr.find_first_of("(") + 1, vecStr.find_last_of(")") - 1);

    vector<string> components;
    split(components, vecStr, is_any_of(", \t"), token_compress_on);

    x = stod(components[0]);
    y = stod(components[1]);
}

Scalar Vector2D::mag() const
{
    return sqrt(x * x + y * y);
}

Scalar Vector2D::magSqr() const
{
    return x * x + y * y;
}

Vector2D Vector2D::abs() const
{
    return Vector2D(std::abs(x), std::abs(y));
}

Vector2D Vector2D::unitVec() const
{
    Scalar invMag = 1. / mag();
    return Vector2D(x * invMag, y * invMag);
}

Vector2D Vector2D::normalVec() const
{
    return Vector2D(y, -x);
}

Vector2D Vector2D::tangentVec() const
{
    return Vector2D(-y, x);
}

Scalar Vector2D::angle() const
{
    return atan2(y, x);
}

Scalar Vector2D::angle(const Vector2D &other) const
{
    return atan2(y - other.y, x - other.x);
}

bool Vector2D::isParallel(const Vector2D &other) const
{
    return fabs(dot(*this, other) * dot(*this, other) - magSqr() * other.magSqr()) < EPSILON_;
}

Vector2D Vector2D::rotate(Scalar theta) const
{
    Scalar ct = std::cos(theta), st = std::sin(theta);
    return Vector2D(x * ct - y * st, x * st + y * ct);
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

Scalar &Vector2D::operator()(int component)
{
    switch (component)
    {
        case 0:
            return x;
        case 1:
            return y;
        default:
            throw Exception("Vector2D", "operator()", "invalid component.");
    }
}

Scalar Vector2D::operator()(int component) const
{
    switch (component)
    {
        case 0:
            return x;
        case 1:
            return y;
        default:
            throw Exception("Vector2D", "operator()", "invalid component.");
    }
}

Vector2D &Vector2D::operator+=(const Vector2D &other)
{
    x += other.x;
    y += other.y;
    return *this;
}

Vector2D &Vector2D::operator-=(const Vector2D &other)
{
    x -= other.x;
    y -= other.y;
    return *this;
}

Vector2D &Vector2D::operator*=(Scalar other)
{
    x *= other;
    y *= other;
    return *this;
}

Vector2D &Vector2D::operator/=(Scalar other)
{
    x /= other;
    y /= other;
    return *this;
}

bool Vector2D::operator<(const Vector2D &rhs) const
{
    return x < rhs.x || (x == rhs.x && y < rhs.y);
}

//- External functions

std::ostream &operator<<(std::ostream &os, const Vector2D &vec)
{
    return os << std::to_string(vec);
}

Vector2D operator+(Vector2D lhs, const Vector2D &rhs)
{
    lhs += rhs;
    return lhs;
}

Vector2D operator-(Vector2D lhs, const Vector2D &rhs)
{
    lhs -= rhs;
    return lhs;
}

Vector2D operator-(const Vector2D &rhs)
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
    return rhs * lhs;
}

Vector2D operator/(Vector2D lhs, Scalar rhs)
{
    lhs /= rhs;
    return lhs;
}

Scalar dot(const Vector2D &u, const Vector2D &v)
{
    return u.x * v.x + u.y * v.y;
}

Scalar cross(const Vector2D &u, const Vector2D &v)
{
    return u.x * v.y - u.y * v.x;
}

Vector2D pointwise(const Vector2D &u, const Vector2D &v)
{
    return Vector2D(u.x * v.x, u.y * v.y);
}

std::string std::to_string(const Vector2D &vec)
{
    return "(" + std::to_string(vec.x) + ", " + std::to_string(vec.y) + ")";
}



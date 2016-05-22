#ifndef VECTOR_2D_H
#define VECTOR_2D_H

#include <ostream>
#include <math.h>

#include "Types.h"

class Vector2D
{
public:

    Vector2D(Scalar x = 0., Scalar y = 0.);
    Vector2D(std::string vecStr);

    Scalar magSqr() const;
    Scalar mag() const;

    Vector2D unitVec() const;
    Vector2D normalVec() const;
    Vector2D tangentVec() const;
    Vector2D abs() const { return Vector2D(fabs(x), fabs(y)); }
    Scalar angle() const;
    Scalar angle(const Vector2D &other) const;
    bool isParallel(const Vector2D &other) const;

    Vector2D rotate(Scalar theta) const;
    Vector2D transform(const Vector2D& uPrime) const;

    std::string toString() const;

    // Operators
    Scalar& operator()(int component);
    Scalar operator ()(int component) const;
    Vector2D& operator +=(const Vector2D& other);
    Vector2D& operator -=(const Vector2D& other);
    Vector2D& operator *=(Scalar other);
    Vector2D& operator /=(Scalar other);
    bool operator==(const Vector2D& rhs) const;

    bool operator<(const Vector2D& rhs);

    Scalar x, y;

private:

    static const Scalar EPSILON_;
};

std::ostream& operator<<(std::ostream& os, const Vector2D& vec);
Vector2D operator+(Vector2D lhs, const Vector2D& rhs);
Vector2D operator-(Vector2D lhs, const Vector2D& rhs);
Vector2D operator-(const Vector2D& rhs);
Vector2D operator*(Vector2D lhs, Scalar rhs);
Vector2D operator*(Scalar lhs, Vector2D rhs);
Vector2D operator/(Vector2D lhs, Scalar rhs);

Scalar dot(const Vector2D& u, const Vector2D& v);
Scalar cross(const Vector2D& u, const Vector2D& v);
std::string to_string(const Vector2D& vec);

#endif

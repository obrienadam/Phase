#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include "Types.h"

class Vector : public std::vector<Scalar>
{
public:
    Vector(size_t size = 0, Scalar val = 0.);
    Scalar& operator()(size_t i) { return std::vector<Scalar>::operator [](i); }
    const Scalar& operator()(size_t i) const { return std::vector<Scalar>::operator [](i); }

    Vector& operator+=(const Vector& rhs);
    Vector& operator-=(const Vector& rhs);
    Vector& operator*=(Scalar rhs);
};

Vector operator+(Vector lhs, const Vector& rhs);
Vector operator-(Vector lhs, const Vector& rhs);
Vector operator*(Scalar lhs, Vector rhs);
Vector operator*(Vector lhs, Scalar rhs);

#endif

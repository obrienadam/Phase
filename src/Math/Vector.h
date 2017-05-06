#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include "Types.h"

class Vector : public std::vector<Scalar>
{
public:
    Vector(size_t size = 0, Scalar val = 0.);

    Vector(const Vector &rhs) = default;

    Vector(Vector &&rhs) = default;

    Scalar &operator()(size_t i)
    { return std::vector<Scalar>::operator[](i); }

    const Scalar &operator()(size_t i) const
    { return std::vector<Scalar>::operator[](i); }

    Vector &operator=(const Vector &rhs) = default;

    Vector &operator=(Vector &&rhs) = default;

    Vector &operator+=(const Vector &rhs);

    Vector &operator-=(const Vector &rhs);

    Vector &operator*=(Scalar rhs);

    Vector operator-() const;

    void zero() { std::fill(begin(), end(), 0.); }
};

Vector operator+(Vector lhs, const Vector &rhs);

Vector operator-(Vector lhs, const Vector &rhs);

Vector operator*(Scalar lhs, Vector rhs);

Vector operator*(Vector lhs, Scalar rhs);

#endif

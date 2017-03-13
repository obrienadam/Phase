#include <algorithm>
#include "Vector.h"

Vector::Vector(size_t size, Scalar val)
        :
        std::vector<Scalar>(size, val)
{

}

Vector &Vector::operator+=(const Vector &rhs)
{
    std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), std::plus<Scalar>());
    return *this;
}

Vector &Vector::operator-=(const Vector &rhs)
{
    std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), std::minus<Scalar>());
    return *this;
}

Vector &Vector::operator*=(Scalar rhs)
{
    std::transform(this->begin(), this->end(), this->begin(), [rhs](Scalar val) { return val * rhs; });
    return *this;
}

//- External functions
Vector operator+(Vector lhs, const Vector &rhs)
{
    return lhs += rhs;
}

Vector operator-(Vector lhs, const Vector &rhs)
{
    return lhs -= rhs;
}

Vector operator*(Scalar lhs, Vector rhs)
{
    return rhs *= lhs;
}

Vector operator*(Vector lhs, Scalar rhs)
{
    return lhs *= rhs;
}

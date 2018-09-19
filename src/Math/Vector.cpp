#include <algorithm>
#include "Vector.h"

Vector::Vector(size_t size, Scalar val)
        :
        data_(size, val)
{

}

Vector &Vector::operator+=(const Vector &rhs)
{
    std::transform(data_.begin(), data_.end(), rhs.data_.begin(), data_.begin(), std::plus<Scalar>());
    return *this;
}

Vector &Vector::operator-=(const Vector &rhs)
{
    std::transform(data_.begin(), data_.end(), rhs.data_.begin(), data_.begin(), std::minus<Scalar>());
    return *this;
}

Vector &Vector::operator+=(Scalar rhs)
{
    std::for_each(data_.begin(), data_.end(), [rhs](Scalar &val) { val += rhs; });
    return *this;
}

Vector &Vector::operator-=(Scalar rhs)
{
    std::for_each(data_.begin(), data_.end(), [rhs](Scalar &val) { val -= rhs; });
    return *this;
}

Vector &Vector::operator*=(Scalar rhs)
{
    std::for_each(data_.begin(), data_.end(), [rhs](Scalar &val) { val *= rhs; });
    return *this;
}

Vector &Vector::operator/=(Scalar rhs)
{
    std::for_each(data_.begin(), data_.end(), [rhs](Scalar &val) { val /= rhs; });
    return *this;
}

Vector Vector::operator-() const
{
    Vector newVector(*this);
    std::for_each(newVector.data_.begin(), newVector.data_.end(), [](Scalar &val) { val = -val; });
    return newVector;
}

void Vector::zero()
{
    std::fill(data_.begin(), data_.end(), 0.);
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

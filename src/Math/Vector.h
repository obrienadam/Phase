#ifndef PHASE_VECTOR_H
#define PHASE_VECTOR_H

#include <vector>

#include "Types/Types.h"

class Vector {
public:
  Vector(Size size = 0, Scalar val = 0.);

  Scalar &operator()(Index i) { return data_[i]; }

  Scalar operator()(Index i) const { return data_[i]; }

  size_t size() const { return data_.size(); }

  void resize(Size size) { data_.resize(size); }

  void resize(Size size, Scalar val) { data_.resize(size, val); }

  void clear() { data_.clear(); }

  const std::vector<Scalar> &data() const { return data_; }

  Vector &operator+=(const Vector &rhs);

  Vector &operator-=(const Vector &rhs);

  Vector &operator+=(Scalar rhs);

  Vector &operator-=(Scalar rhs);

  Vector &operator*=(Scalar rhs);

  Vector &operator/=(Scalar rhs);

  Vector operator-() const;

  void zero();

private:
  std::vector<Scalar> data_;
};

Vector operator+(Vector lhs, const Vector &rhs);

Vector operator-(Vector lhs, const Vector &rhs);

Vector operator*(Scalar lhs, Vector rhs);

Vector operator*(Vector lhs, Scalar rhs);

#endif

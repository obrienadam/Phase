#ifndef PHASE_POLY_1D_H
#define PHASE_POLY_1D_H

#include "StaticMatrix.h"
#include "Types/Types.h"

template <Size N> class Poly1D {
public:
  Poly1D(const std::initializer_list<Scalar> &x,
         const std::initializer_list<Scalar> &rhs);

  Scalar operator()(Scalar x) const;

private:
  StaticMatrix _coeffs;
};

#include "Poly1D.tpp"

#endif

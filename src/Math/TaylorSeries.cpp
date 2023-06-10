#include <math.h>

#include "Factorial.h"
#include "Matrix.h"
#include "TaylorSeries.h"

void TaylorSeries::computeTaylorCoeffs(const std::vector<Scalar> &offsets,
                                       int derivOrder) {
  _offsets = offsets;

  Matrix A(_offsets.size(), _offsets.size());

  for (int j = 0; j < _offsets.size(); ++j)
    for (int i = 0; i < _offsets.size(); ++i)
      A(j, i) = std::pow(_offsets[i], j) / factorial(j);

  Matrix b(_offsets.size(), 1);

  for (int i = 0; i < _offsets.size(); ++i)
    b(i, 0) = i == derivOrder ? 1. : 0.;

  A.solve(b);
  _coeffs.assign(b.begin(), b.end());
}

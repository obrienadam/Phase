#ifndef PHASE_TAYLOR_SERIES_H
#define PHASE_TAYLOR_SERIES_H

#include <vector>

#include "Types/Types.h"

class TaylorSeries {
public:
  TaylorSeries() {}

  TaylorSeries(const std::vector<Scalar> &offsets, int derivOrder) {
    computeTaylorCoeffs(offsets, derivOrder);
  }

  void computeTaylorCoeffs(const std::vector<Scalar> &offsets, int derivOrder);

  Size size() const { return _coeffs.size(); }

  int derivOrder() const { return _derivOrder; }

  const std::vector<Scalar> &offsets() const { return _offsets; }

  const std::vector<Scalar> &coeffs() const { return _coeffs; }

protected:
  int _derivOrder;

  std::vector<Scalar> _offsets, _coeffs;
};

#endif

#ifndef PHASE_FACE_STENCIL_H
#define PHASE_FACE_STENCIL_H

#include "Geometry/Point3D.h"

#include "Math/TaylorSeries.h"

#include "Face.h"

class Cell;

class FaceStencil {
public:
  FaceStencil(const Cell &cell, Face::Direction dir, int forwardBias,
              int backwardBias);

  FaceStencil(const Cell &cell, Face::Direction dir, int order);

  Size maxForwardShift();

  Size maxBackwardShift();

  const Cell &operator()(Size shift) const;

  const std::vector<Ref<const Cell>> &cells() const { return _cells; }

  const Face &face() const { return _face; }

  const std::vector<Scalar> &coeffs() const { return _taylorCoeffs.coeffs(); }

protected:
  Point3D _xf;

  Vector3D _sf, _rf;

  const Cell &_cell;

  const Face &_face;

  Face::Direction _dir;

  std::vector<Ref<const Cell>> _cells;

  TaylorSeries _taylorCoeffs;
};

#endif

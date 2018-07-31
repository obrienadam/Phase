#include "BoundaryFaceStencil.h"

BoundaryFaceStencil::BoundaryFaceStencil(const Cell &cell, Face::Direction dir, int order, BoundaryType btype)
    :
      FaceStencil(cell, dir, order),
      _btype(btype)
{
    _faceCoeffs = _taylorCoeffs.coeffs();

    for(Scalar &c: _faceCoeffs)
        c /= _taylorCoeffs.coeffs().back();

    _faceCoeffs.pop_back();
}

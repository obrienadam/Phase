#include "BoundaryFaceStencil.h"

BoundaryFaceStencil::BoundaryFaceStencil(const Cell &cell, Face::Direction dir, int order)
    :
      FaceStencil(cell, dir, order)
{
    _faceCoeffs = _taylorCoeffs.coeffs();

    for(Scalar &c: _faceCoeffs)
        c /= _taylorCoeffs.coeffs().back();

    _faceCoeffs.pop_back();
}

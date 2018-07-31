#ifndef PHASE_BOUNDARY_FACE_STENCIL_H
#define PHASE_BOUNDARY_FACE_STENCIL_H

#include "FaceStencil.h"

class BoundaryFaceStencil: public FaceStencil
{
public:

    enum BoundaryType{FIXED, ZERO_GRADIENT};

    BoundaryFaceStencil(const Cell& cell, Face::Direction dir, int order, BoundaryType btype);

    BoundaryType btype() const
    { return _btype; }

    const std::vector<Scalar> &faceCoeffs() const
    { return _faceCoeffs; }

protected:

    BoundaryType _btype;

    std::vector<Scalar> _faceCoeffs;
};

#endif

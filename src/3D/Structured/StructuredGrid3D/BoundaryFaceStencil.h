#ifndef PHASE_BOUNDARY_FACE_STENCIL_H
#define PHASE_BOUNDARY_FACE_STENCIL_H

#include "FaceStencil.h"

class BoundaryFaceStencil: public FaceStencil
{
public:

    BoundaryFaceStencil(const Cell& cell, Face::Direction dir, int order);

    const std::vector<Scalar> &faceCoeffs() const
    { return _faceCoeffs; }

protected:

    std::vector<Scalar> _faceCoeffs;
};

#endif

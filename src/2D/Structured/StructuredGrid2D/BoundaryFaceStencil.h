#ifndef PHASE_BOUNDARY_FACE_STENCIL_H
#define PHASE_BOUNDARY_FACE_STENCIL_H

#include "Stencil.h"

class Face;

class BoundaryFaceStencil: public Stencil
{
public:

    BoundaryFaceStencil(const Cell &cell, Coordinates::Direction zeta);

    const Vector2D& rf() const
    { return _rf; }

    const Vector2D& sf() const
    { return _sf; }

protected:

    const Face &_face;

    Vector2D _rf, _sf;

};

#endif

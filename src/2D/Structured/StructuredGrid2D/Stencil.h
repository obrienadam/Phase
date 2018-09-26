#ifndef PHASE_STENCIL_H
#define PHASE_STENCIL_H

#include "Geometry/Vector2D.h"

#include "Coordinates.h"

class Cell;

class Stencil
{
public:

    Stencil(const Cell &cell, Coordinates::Direction zeta);

    const Cell &cell(int inc) const;

    Vector2D rc(int inc) const;

    Coordinates::Direction zeta() const
    { return _zeta; }

protected:

    const Cell &_cell;

    Coordinates::Direction _zeta;
};

#endif

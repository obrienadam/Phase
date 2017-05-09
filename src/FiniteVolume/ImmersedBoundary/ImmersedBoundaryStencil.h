#ifndef IMMERSED_BOUNDARY_STENCIL_H
#define IMMERSED_BOUNDARY_STENCIL_H

#include "Cell.h"

class ImmersedBoundaryStencil
{
public:

    ImmersedBoundaryStencil(const Cell& cell) : cell_(cell) {}

protected:

    const Cell &cell_;
};

#endif

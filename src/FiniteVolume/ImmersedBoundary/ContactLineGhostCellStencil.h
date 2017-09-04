#ifndef CONTACT_LINE_GHOST_CELL_STENCIL_H
#define CONTACT_LINE_GHOST_CELL_STENCIL_H

#include "GhostCellStencil.h"

class ContactLineGhostCellStencil: public GhostCellStencil
{
public:
    ContactLineGhostCellStencil(const Cell &cell,
                                const GhostCellImmersedBoundaryObject &ibObj,
                                const ScalarFiniteVolumeField& gamma,
                                Scalar theta);
};


#endif

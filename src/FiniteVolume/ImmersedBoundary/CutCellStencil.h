#ifndef CUT_CELL_STENCIL_H
#define CUT_CELL_STENCIL_H

#include "Cell.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

class CutCellStencil
{
public:

    CutCellStencil(const Cell& cell, const Shape2D& shape);

private:

    const Cell& cell_;

};

#endif

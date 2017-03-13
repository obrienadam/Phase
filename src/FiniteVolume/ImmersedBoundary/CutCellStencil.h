#ifndef CUT_CELL_STENCIL_H
#define CUT_CELL_STENCIL_H

#include "Cell.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

class CutCellStencil
{
public:

    CutCellStencil(const Cell &cell, const Shape2D &shape);

    const std::vector<LineSegment2D> &faces() const
    { return faces_; }

    const std::vector<LineSegment2D> &ibFaces() const
    { return ibFaces_; }

private:

    const Cell &cell_;
    std::vector< Ref<const Cell> > nbCells_;
    std::vector<LineSegment2D> faces_, ibFaces_;

};

#endif

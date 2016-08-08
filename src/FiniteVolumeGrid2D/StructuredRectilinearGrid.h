#ifndef STRUCTURED_RECTILINEAR_GRID_H
#define STRUCTURED_RECTILINEAR_GRID_H

#include "FiniteVolumeGrid2D.h"
#include "Input.h"

class StructuredRectilinearGrid : public FiniteVolumeGrid2D
{
public:

    StructuredRectilinearGrid(Scalar width, Scalar height, Size nCellsX, Size nCellsY);

    Cell& operator ()(Label i, Label j);
    const Cell& operator()(Label i, Label j) const;
    const Node& node(Label i, Label j) const;

protected:

    Size nCellsX_, nCellsY_;
    Scalar width_, height_;

};

#endif

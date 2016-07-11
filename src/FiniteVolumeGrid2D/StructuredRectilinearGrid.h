#ifndef STRUCTURED_RECTILINEAR_GRID_H
#define STRUCTURED_RECTILINEAR_GRID_H

#include "FiniteVolumeGrid2D.h"
#include "Input.h"

class StructuredRectilinearGrid : public FiniteVolumeGrid2D
{
public:

    StructuredRectilinearGrid(Scalar width, Scalar height, int nCellsX, int nCellsY);

    Cell& operator ()(int i, int j);
    const Cell& operator()(int i, int j) const;
    const Node& node(int i, int j) const;

protected:

    int nCellsX_, nCellsY_;
    Scalar width_, height_;

};

#endif

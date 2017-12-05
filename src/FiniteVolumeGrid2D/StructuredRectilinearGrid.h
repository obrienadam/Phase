#ifndef STRUCTURED_RECTILINEAR_GRID_H
#define STRUCTURED_RECTILINEAR_GRID_H

#include <set>

#include "FiniteVolumeGrid2D.h"
#include "Input.h"

class StructuredRectilinearGrid : public FiniteVolumeGrid2D
{
public:

    StructuredRectilinearGrid(Scalar width,
                              Scalar height,
                              Size nCellsX,
                              Size nCellsY,
                              Scalar convertToMeters,
                              const std::vector<std::pair<Scalar, Scalar>> &xDimRefinements,
                              const std::vector<std::pair<Scalar, Scalar>> &yDimRefinements,
                              const Point2D& origin);

    Cell &operator()(Label i, Label j);

    const Cell &operator()(Label i, Label j) const;

    const Node &node(Label i, Label j) const;

protected:

    void refineDims(Scalar start, Scalar end, std::vector<Scalar> &dims);

    Size nCellsX_, nCellsY_;
    Scalar width_, height_;

};

#endif

#ifndef STRUCTURED_RECTILINEAR_GRID_H
#define STRUCTURED_RECTILINEAR_GRID_H

#include "System/Input.h"

#include "FiniteVolumeGrid2D.h"

class StructuredRectilinearGrid : public FiniteVolumeGrid2D
{
public:

    StructuredRectilinearGrid(const Input& input);

    void init(Scalar width,
              Scalar height,
              Size nCellsX,
              Size nCellsY,
              Scalar convertToMeters,
              const std::vector<Point2D> &xDimRefinements,
              const std::vector<Point2D> &yDimRefinements,
              const Point2D &origin);

    Cell &operator()(Label i, Label j);

    const Cell &operator()(Label i, Label j) const;

    const Node &node(Label i, Label j) const;

    bool isEquidistant() const
    { return std::abs(width_ / nCellsX_ - height_ / nCellsY_) < 1e-12; }

    Scalar hx() const
    { return width_ / nCellsX_; }

    Scalar hy() const
    { return height_ / nCellsY_; }

    Scalar h() const;

protected:

    std::vector<Scalar> refineDims(Scalar start, Scalar end, const std::vector<Scalar> &dims);

    Size nCellsX_, nCellsY_;

    Scalar width_, height_;

};

#endif

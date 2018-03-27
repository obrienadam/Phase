#ifndef PHASE_STRUCTURED_GRID_2D_H
#define PHASE_STRUCTURED_GRID_2D_H

#include <vector>

#include "Geometry/Point2D.h"

class StructuredGrid2D
{
public:

    class Index
    {
    public:
        size_t i, j, id;
    };

    StructuredGrid2D(size_t nNodesI, size_t nNodesJ, Scalar lx, Scalar ly);

    size_t id(size_t i, size_t j) const
    { return j * nNodesI_ + i; }

    size_t nNodes() const
    { return nNodesI_ * nNodesJ_; }

    const Point2D &node(size_t i, size_t j) const;

    const Point2D &node(size_t id) const;

    Scalar dxe(size_t i, size_t j) const;

    Scalar dxw(size_t i, size_t j) const;

    Scalar dxn(size_t i, size_t j) const;

    Scalar dxs(size_t i, size_t j) const;

protected:

    size_t nNodesI_, nNodesJ_;

    Scalar lx_, ly_;

    std::vector<Point2D> nodes_;
};


#endif

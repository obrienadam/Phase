#ifndef HIGH_ORDER_STENCIL_H
#define HIGH_ORDER_STENCIL_H

#include "ImmersedBoundary.h"

class HighOrderStencil
{
public:
    HighOrderStencil(const Cell &cell, const ImmersedBoundary &ib);

    HighOrderStencil(const Cell &cell, const ImmersedBoundaryObject &ibObj);

    const std::vector<Ref<const Cell>> &cells() const
    { return cells_; }

    const std::vector<Ref<const Cell>> &compatCells() const
    { return compatCells_; }

    const std::vector<Point2D> &compatPoints() const
    { return compatPoints_; }

    const std::vector<Point2D> &boundaryPoints() const
    { return boundaryPoints_; }

    size_t numEqns() const
    { return cells_.size() + 2 * compatCells().size() + boundaryPoints_.size(); }

private:

    std::vector<Ref<const Cell>> cells_, compatCells_;
    std::vector<Point2D> compatPoints_;
    std::vector<Point2D> boundaryPoints_;
};

#endif
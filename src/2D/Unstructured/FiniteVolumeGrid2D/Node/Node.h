#ifndef PHASE_NODE_H
#define PHASE_NODE_H

#include "Geometry/Point2D.h"

class Cell;

class Face;

class FiniteVolumeGrid2D;

class Node : public Point2D
{
public:

    explicit Node(Scalar x, Scalar y, const FiniteVolumeGrid2D &grid);

    explicit Node(const Point2D &point, const FiniteVolumeGrid2D &grid);

    Label id() const
    { return id_; }

    void setId(Label id)
    { id_ = id; }

    void addCell(const Cell &cell)
    { cells_.emplace_back(cell); }

    const std::vector<Ref<const Cell>> &cells() const
    { return cells_; }

    bool isBoundaryNode() const;

    std::vector<Scalar> volumeWeights() const;

    std::vector<Scalar> distanceWeights() const;

    const Point2D &centroid() const //- For compatibility with Groups
    { return *this; }

    bool cellBounded(const Point2D& pt) const;

protected:

    Label id_;

    std::vector<Ref<const Cell>> cells_;

    const FiniteVolumeGrid2D &grid_;
};

#endif

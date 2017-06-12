#ifndef NODE_H
#define NODE_H

#include "Point2D.h"

class Cell;

class Face;

class FiniteVolumeGrid2D;

class Node : public Point2D
{
public:

    Node(Scalar x, Scalar y, const FiniteVolumeGrid2D &grid);

    explicit Node(const Point2D &point, const FiniteVolumeGrid2D &grid);

    Label id() const
    { return id_; }

    void setId(Label id)
    { id_ = id; }

    void addCell(const Cell &cell);

    const std::vector<Label> &cellIds() const
    { return cellIds_; }

    const std::vector<Ref<const Cell>> cells() const;

    const Point2D &centroid() const //- For compatibility with Groups
    { return *this; }

protected:

    Label id_;

    std::vector<Label> cellIds_;
    const std::vector<Cell> &cells_;
};

#endif

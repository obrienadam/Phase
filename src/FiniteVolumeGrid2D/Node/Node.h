#ifndef NODE_H
#define NODE_H

#include "Point2D.h"

class Cell;
class Face;

class Node : public Point2D
{
public:

    Node(Scalar x, Scalar y, Label id) : Point2D(x, y), id_(id) { }
    explicit Node(const Point2D& point, Label id) : Point2D(point), id_(id) {}

    Label id() const { return id_; }

    const std::vector<Ref<const Cell>>& cells() const { return cells_; }
    const std::vector<Ref<const Face>>& faces() const { return faces_; }

protected:

    Label id_;

    std::vector<Ref<const Cell>> cells_;
    std::vector<Ref<const Face>> faces_;

    friend class FiniteVolumeGrid2D;
};

#endif

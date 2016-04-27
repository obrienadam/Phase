#ifndef NODE_H
#define NODE_H

#include "Point2D.h"

class Node : public Point2D
{
public:

    Node(Scalar x = 0., Scalar y = 0., size_t id = 0) : Point2D(x, y), id_(id) { }

    size_t id() const { return id_; }

    bool operator ==(const Node& other) const;

private:

    size_t id_;

    friend class FiniteVolumeGrid2D;
};

#endif

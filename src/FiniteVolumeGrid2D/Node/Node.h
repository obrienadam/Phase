#ifndef NODE_H
#define NODE_H

#include "Point2D.h"

class Node : public Point2D
{
public:

    Node(Scalar x, Scalar y, Label id) : Point2D(x, y), id_(id) { }
    explicit Node(const Point2D& point, Label id) : Point2D(point), id_(id) {}

    Label id() const { return id_; }

private:

    Label id_;
};

#endif

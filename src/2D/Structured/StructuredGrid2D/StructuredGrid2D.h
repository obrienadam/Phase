#ifndef PHASE_STRUCTURED_GRID_2D_H
#define PHASE_STRUCTURED_GRID_2D_H

#include <vector>

#include "System/Communicator.h"

#include "Geometry/Point2D.h"

class StructuredGrid2D
{
public:

    enum Boundary
    {
        EAST, WEST, NORTH, SOUTH
    };

    class Node : public Point2D
    {
    public:

        Node(const Point2D &loc, Label i, Label j, Label id)
                : Point2D(loc), _i(i), _j(j), _id(id)
        {}

        Label i() const
        { return _i; }

        Label j() const
        { return _j; }

        Label id() const
        { return _id; }

    protected:

        Label _i, _j, _id;
    };

    StructuredGrid2D()
    {}

    StructuredGrid2D(Size nNodesI, Size nNodesJ, Scalar lx, Scalar ly);

    void init(Size nNodesI, Size nNodesJ, Scalar lx, Scalar ly);

    size_t id(size_t i, size_t j) const
    { return j * nNodesI_ + i; }

    size_t nNodesI() const
    { return nNodesI_; }

    size_t nNodesJ() const
    { return nNodesJ_; }

    size_t nNodes() const
    { return nNodesI_ * nNodesJ_; }

    const Point2D &node(size_t i, size_t j) const
    { return nodes_[j * nNodesI_ + i]; }

    Scalar dxe(size_t i, size_t j) const;

    Scalar dxw(size_t i, size_t j) const;

    Scalar dxn(size_t i, size_t j) const;

    Scalar dxs(size_t i, size_t j) const;

    const Communicator &comm() const
    { return _comm; }

protected:

    size_t nNodesI_, nNodesJ_;

    Scalar lx_, ly_;

    std::vector<Node> nodes_;

    std::vector<int> ownership_, globalIds_;

    Communicator _comm;
};


#endif

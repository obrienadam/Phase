#ifndef PHASE_STRUCTURED_GRID_2D_H
#define PHASE_STRUCTURED_GRID_2D_H

#include <vector>

#include "System/Communicator.h"

#include "Geometry/Point2D.h"

class StructuredGrid2D
{
public:

    enum Boundary{EAST, WEST, NORTH, SOUTH};

    class IndexSet
    {
    public:
        size_t ibegin, iend, jbegin, jend;
    };

    StructuredGrid2D(size_t nNodesI, size_t nNodesJ, Scalar lx, Scalar ly);

    size_t id(size_t i, size_t j) const
    { return j * nNodesI_ + i; }

    size_t nNodesI() const
    { return nNodesI_; }

    size_t nNodesJ() const
    { return nNodesJ_; }

    size_t nCellsI() const
    { return nCellsI_; }

    size_t nCellsJ() const
    { return nCellsJ_; }

    size_t iend() const
    { return nCellsI_ - 1; }

    size_t jend() const
    { return nCellsJ_ - 1; }

    size_t nNodes() const
    { return nNodesI_ * nNodesJ_; }

    const Point2D &node(size_t i, size_t j) const
    { return nodes_[j * nNodesI_ + i]; }

    Scalar vol(size_t i, size_t j) const;

    Scalar dxe(size_t i, size_t j) const;

    Scalar dxw(size_t i, size_t j) const;

    Scalar dxn(size_t i, size_t j) const;

    Scalar dxs(size_t i, size_t j) const;

    const Communicator& comm() const
    { return _comm; }

protected:

    size_t nNodesI_, nNodesJ_;

    size_t nCellsI_, nCellsJ_;

    size_t nIFacesI_, nIFacesJ_;

    size_t nJFacesI_, nJFacesJ_;

    Scalar lx_, ly_;

    std::vector<Point2D> nodes_;

    std::vector<int> ownership_, globalIds_;

    IndexSet localCells_;

    Communicator _comm;
};


#endif

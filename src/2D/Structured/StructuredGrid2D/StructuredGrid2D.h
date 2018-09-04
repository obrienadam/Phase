#ifndef PHASE_STRUCTURED_GRID_2D_H
#define PHASE_STRUCTURED_GRID_2D_H

#include "System/Input.h"
#include "System/Communicator.h"

#include "Cell.h"
#include "Face.h"
#include "Set.h"

class StructuredGrid2D
{
public:

    StructuredGrid2D() : _comm(std::make_shared<Communicator>())
    {}

    StructuredGrid2D(Size nCellsI, Size nCellsJ, Scalar lx, Scalar ly);

    StructuredGrid2D(const std::vector<Scalar> &xcoords, const std::vector<Scalar> &ycoords);

    StructuredGrid2D(const Input &input);

    void init(Size nCellsI, Size nCellsJ, Scalar lx, Scalar ly);

    void init(std::vector<Scalar> xcoords, std::vector<Scalar> ycoords);

    //- Parameters
    Size nNodesI() const
    { return _nCellsI + 1; }

    Size nNodesJ() const
    { return _nCellsJ + 1; }

    Size nCellsI() const
    { return _nCellsI; }

    Size nCellsJ() const
    { return _nCellsJ; }

    Size nFacesI() const
    { return nNodesI() * nCellsJ(); }

    Size nFacesJ() const
    { return nCellsI() * nNodesJ(); }

    //- Element access
    const std::vector<Cell> &cells() const
    { return _cells; }

    const Cell& operator()(Label i, Label j) const
    { return _cells[j * _nCellsI + i]; }

    const Cell &cell(const Cell& cell, Orientation dir, int offset) const;

    const Set<Cell> &localCells() const
    { return _localCells; }

    const std::vector<Point2D> nodes() const
    { return _nodes; }

    const Point2D &node(Label i, Label j) const
    { return _nodes[j * nNodesI() + i]; }

    const std::vector<Face> &faces() const
    { return _faces; }

    const std::vector<Face> &ifaces() const
    { return _ifaces; }

    const std::vector<Face> &jfaces() const
    { return _jfaces; }

    const Face &iface(Label i, Label j) const
    { return _ifaces[j * nNodesI() + i]; }

    const Face &jface(Label i, Label j) const
    { return _jfaces[i * nNodesJ() + j]; }

    const Face &face(const Cell &cell, Orientation dir) const;

    //- Misc mesh functions

    Scalar dh(const Cell &cell, Orientation dir, int offset) const;

    //- Communication
    const Communicator &comm() const
    { return *_comm; }

    const std::vector<int> &ownership() const
    { return _ownership; }

protected:

    //- Mesh parameters
    Size _nCellsI, _nCellsJ;

    Scalar _lx, _ly;

    //- Mesh entities
    std::vector<Point2D> _nodes;

    std::vector<Cell> _cells;

    std::vector<Face> _faces, _ifaces, _jfaces;

    //- Parallel
    std::shared_ptr<const Communicator> _comm;

    Set<Cell> _localCells;

    std::vector<int> _ownership;
};


#endif

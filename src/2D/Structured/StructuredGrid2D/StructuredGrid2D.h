#ifndef PHASE_STRUCTURED_GRID_2D_H
#define PHASE_STRUCTURED_GRID_2D_H

#include <vector>

#include "System/Communicator.h"

#include "Cell.h"
#include "Face.h"

class StructuredGrid2D
{
public:

    enum Coordinate{I, J};

    enum CoordinateDirection{I_POS, I_NEG, J_POS, J_NEG};

    StructuredGrid2D()
    {}

    StructuredGrid2D(Size nCellsI, Size nCellsJ, Scalar lx, Scalar ly);

    void init(Size nCellsI, Size nCellsJ, Scalar lx, Scalar ly);

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

    const Cell &cell(Label i, Label j) const
    { return _cells[j * _nCellsI + i]; }

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

    const Communicator &comm() const
    { return *_comm; }

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

    std::vector<int> _ownership;
};


#endif

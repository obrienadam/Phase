#ifndef PHASE_STRUCTURED_GRID_3D_H
#define PHASE_STRUCTURED_GRID_3D_H

#include <memory>

#include "System/Input.h"
#include "System/Communicator.h"
#include "CellSet.h"
#include "FaceSet.h"

class StructuredGrid3D
{
public:

    enum Index{I, J, K};

    StructuredGrid3D(Size nCellsI = 0, Size nCellsJ = 1, Size nCellsK = 1, Scalar lx = 1., Scalar ly = 1., Scalar lz = 1.);

    StructuredGrid3D(const Input &input);

    //- Grid info
    Scalar lx() const
    { return _lx; }

    Scalar ly() const
    { return _ly; }

    Scalar lz() const
    { return _lz; }

    Size nNodesI() const
    { return _nCellsI + 1; }

    Size nNodesJ() const
    { return _nCellsJ + 1; }

    Size nNodesK() const
    { return _nCellsK + 1; }

    Size nCellsI() const
    { return _nCellsI; }

    Size nCellsJ() const
    { return _nCellsJ; }

    Size nCellsK() const
    { return _nCellsK; }

    Size nCells() const
    { return _nCells; }

    Size iMax() const
    { return _nCellsI - 1; }

    Size jMax() const
    { return _nCellsJ - 1; }

    Size kMax() const
    { return _nCellsK - 1; }

    //- Access methods
    const Point3D &node(Label i, Label j, Label k) const
    { return _nodes[(_nCellsI + 1) * (_nCellsJ + 1) * k + (_nCellsI + 1) * j + i]; }

    const Face &operator()(Label i, Label j, Label k, Face::Direction f) const;

    const Face &operator()(const Cell &cell, Face::Direction f) const
    { return operator()(cell.i(), cell.j(), cell.k(), f); }

    const std::vector<Point3D> &nodes() const
    { return _nodes; }

    const std::vector<Cell> &cells() const
    { return _cells; }

    const std::vector<Face> &ifaces() const
    { return _ifaces; }

    const std::vector<Face> &jfaces() const
    { return _jfaces; }

    const std::vector<Face> &kfaces() const
    { return _kfaces; }

    const Cell& operator()(Label i, Label j, Label k) const;

    const Communicator& comm() const
    { return *_comm; }

protected:

    Scalar _lx, _ly, _lz;

    Size _nCellsI, _nCellsJ, _nCellsK, _nCells;

    std::vector<Point3D> _nodes;

    std::vector<Cell> _cells;

    CellSet _localCells;

    std::vector<Face> _ifaces, _jfaces, _kfaces;

    std::shared_ptr<const Communicator> _comm;
};

#endif

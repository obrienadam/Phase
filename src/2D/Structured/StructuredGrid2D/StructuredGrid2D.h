#ifndef PHASE_STRUCTURED_GRID_2D_H
#define PHASE_STRUCTURED_GRID_2D_H

#include <unordered_map>

#include "System/Communicator.h"
#include "System/Input.h"

#include "Cell.h"
#include "Face.h"
#include "Set.h"

class StructuredGrid2D {
public:
  StructuredGrid2D();

  StructuredGrid2D(const StructuredGrid2D &) = delete;

  StructuredGrid2D(Size nCellsI, Size nCellsJ, Scalar lx, Scalar ly,
                   int nBufferCells);

  StructuredGrid2D(const Input &input);

  void init(Size nCellsI, Size nCellsJ, std::pair<Scalar, Scalar> xb,
            std::pair<Scalar, Scalar> yb);

  void init(Size nCellsI, Size nCellsJ, Scalar lx, Scalar ly, int nbuff);

  //- Parameters
  Size nNodesI() const { return _nCellsI + 1; }

  Size nNodesJ() const { return _nCellsJ + 1; }

  Size nCellsI() const { return _nCellsI; }

  Size nCellsJ() const { return _nCellsJ; }

  Size nFacesI() const { return nNodesI() * nCellsJ(); }

  Size nFacesJ() const { return nCellsI() * nNodesJ(); }

  //- Element access
  const std::vector<Cell> &cells() const { return _cells; }

  const Cell &operator()(Label i, Label j) const {
    return _cells[j * _nCellsI + i];
  }

  const Cell &cell(const Cell &cell, Coordinates::Direction dir,
                   int offset) const;

  const Set<Cell> &localCells() const { return _localCells; }

  const std::vector<Point2D> nodes() const { return _nodes; }

  const Point2D &node(Label i, Label j) const {
    return _nodes[j * nNodesI() + i];
  }

  const std::vector<Face> &faces() const { return _faces; }

  const std::vector<Face> &ifaces() const { return _ifaces; }

  const std::vector<Face> &jfaces() const { return _jfaces; }

  const Face &iface(Label i, Label j) const {
    return _ifaces[j * nNodesI() + i];
  }

  const Face &jface(Label i, Label j) const {
    return _jfaces[i * nNodesJ() + j];
  }

  const Face &face(const Cell &cell, Coordinates::Direction dir) const;

  //- Misc mesh functions

  Scalar dh(const Cell &cell, Coordinates::Direction dir, int offset) const;

  Size maxInc(const Cell &cell, Coordinates::Direction dir) const;

  //- Communication and parallel
  const Communicator &comm() const { return *_comm; }

  const std::vector<int> &ownership() const { return _ownership; }

  const Set<Cell> &sendBuffer(int proc) const { return _sendBuffers[proc]; }

  const Set<Cell> &recvBuffer(int proc) const { return _recvBuffers[proc]; }

protected:
  std::array<Size, 2> computeBlockDims(bool stripPartitioning) const;

  void initParallel(const std::vector<int> &ownership,
                    const std::vector<Label> &gids);

  //- Mesh parameters
  Size _nCellsI, _nCellsJ;

  std::pair<Scalar, Scalar> _xb, _yb;

  //- Mesh entities
  std::vector<Point2D> _nodes;

  std::vector<Cell> _cells;

  std::vector<Face> _faces, _ifaces, _jfaces;

  //- Parallel
  std::shared_ptr<const Communicator> _comm;

  std::vector<int> _ownership;

  std::unordered_map<Label, Label> _globalToLocalIdMap;

  Set<Cell> _localCells, _globalCells;

  std::vector<Set<Cell>> _sendBuffers, _recvBuffers;
};

#endif

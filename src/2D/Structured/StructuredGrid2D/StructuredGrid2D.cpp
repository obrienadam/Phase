#include "System/Exception.h"

#include "StructuredGrid2D.h"

StructuredGrid2D::StructuredGrid2D() : _comm(std::make_shared<Communicator>()) {
  _sendBuffers.resize(_comm->nProcs());
  _recvBuffers.resize(_comm->nProcs());
}

StructuredGrid2D::StructuredGrid2D(Size nCellsI, Size nCellsJ, Scalar lx,
                                   Scalar ly, int nBufferCells)
    : StructuredGrid2D() {
  init(nCellsI, nCellsJ, lx, ly, nBufferCells);
}

StructuredGrid2D::StructuredGrid2D(const Input &input)
    : StructuredGrid2D(input.caseInput().get<Size>("Grid.nCellsX"),
                       input.caseInput().get<Size>("Grid.nCellsY"),
                       input.caseInput().get<Scalar>("Grid.width"),
                       input.caseInput().get<Scalar>("Grid.height"),
                       input.caseInput().get<int>("Grid.nBufferCells", 2)) {}

void StructuredGrid2D::init(Size nCellsI, Size nCellsJ,
                            std::pair<Scalar, Scalar> xb,
                            std::pair<Scalar, Scalar> yb) {
  _nCellsI = nCellsI;
  _nCellsJ = nCellsJ;
  _xb = xb;
  _yb = yb;

  Scalar dx = (_xb.second - _xb.first) / _nCellsI;
  Scalar dy = (_yb.second - _yb.first) / _nCellsJ;

  _nodes.clear();

  for (auto j = 0; j < nNodesJ(); ++j)
    for (auto i = 0; i < nNodesI(); ++i)
      _nodes.push_back(Point2D(i * dx + _xb.first, j * dy + _yb.first));

  _cells.clear();

  for (auto j = 0; j < _nCellsJ; ++j)
    for (auto i = 0; i < _nCellsI; ++i)
      _cells.push_back(Cell(*this, i, j));

  _ifaces.clear();

  for (int j = 0; j < _nCellsJ; ++j)
    for (int i = 0; i < nNodesI(); ++i) {
      _ifaces.push_back(Face(*this, Coordinates::I, i, j));
    }

  _jfaces.clear();

  for (int i = 0; i < _nCellsI; ++i)
    for (int j = 0; j < nNodesJ(); ++j) {
      _jfaces.push_back(Face(*this, Coordinates::J, i, j));
    }

  _faces.clear();
  for (const Face &f : _ifaces)
    _faces.push_back(f);

  for (const Face &f : _jfaces)
    _faces.push_back(f);

  _localCells.add(_cells.begin(), _cells.end());
  _ownership.resize(_cells.size(), _comm->rank());
}

void StructuredGrid2D::init(Size nCellsI, Size nCellsJ, Scalar lx, Scalar ly,
                            int nbuff) {
  auto dim = computeBlockDims(false);

  _comm->printf("Grid partitioning dimensions (%d,%d).\n", (int)dim[0],
                (int)dim[1]);

  std::vector<int> ownership(nCellsI * nCellsJ);
  std::array<Size, 2> irange = {0, 0}, jrange = {0, 0};
  std::array<Size, 2> irangeLocal = {0, 0}, jrangeLocal = {0, 0};

  Size qi = nCellsI / dim[0];
  Size ri = nCellsI % dim[0];
  Size qj = nCellsJ / dim[1];
  Size rj = nCellsJ % dim[1];

  for (auto blockJ = 0, proc = 0; blockJ < dim[1]; ++blockJ)
    for (auto blockI = 0; blockI < dim[0]; ++blockI, ++proc) {
      irange[0] = blockI * qi + std::min(ri, (Size)blockI);
      irange[1] = irange[0] + qi + (ri > blockI ? 1 : 0);
      jrange[0] = blockJ * qj + std::min(rj, (Size)blockJ);
      jrange[1] = jrange[0] + qj + (rj > blockJ ? 1 : 0);

      for (auto j = jrange[0]; j < jrange[1]; ++j)
        for (auto i = irange[0]; i < irange[1]; ++i)
          ownership[j * nCellsI + i] = proc;

      if (proc == _comm->rank()) {
        irangeLocal[0] = std::max((int)irange[0] - nbuff, 0);
        irangeLocal[1] = std::min((int)irange[1] + nbuff, (int)nCellsI);
        jrangeLocal[0] = std::max((int)jrange[0] - nbuff, 0);
        jrangeLocal[1] = std::min((int)jrange[1] + nbuff, (int)nCellsJ);
      }
    }

  std::vector<int> localOwnership;
  std::vector<Label> localGlobalIds;

  for (auto j = jrangeLocal[0]; j < jrangeLocal[1]; ++j)
    for (auto i = irangeLocal[0]; i < irangeLocal[1]; ++i) {
      localGlobalIds.emplace_back(j * nCellsI + i);
      localOwnership.emplace_back(ownership[localGlobalIds.back()]);
    }

  Scalar dx = lx / nCellsI;
  Scalar dy = ly / nCellsJ;

  auto xb = std::make_pair(dx * irangeLocal[0], dx * irangeLocal[1]);
  auto yb = std::make_pair(dy * jrangeLocal[0], dy * jrangeLocal[1]);

  init(irangeLocal[1] - irangeLocal[0], jrangeLocal[1] - jrangeLocal[0], xb,
       yb);
  initParallel(localOwnership, localGlobalIds);
}

const Cell &StructuredGrid2D::cell(const Cell &cell, Coordinates::Direction dir,
                                   int offset) const {
  switch (dir) {
  case Coordinates::I_POS:
    return operator()(cell.i() + offset, cell.j());
  case Coordinates::I_NEG:
    return operator()(cell.i() - offset, cell.j());
  case Coordinates::J_POS:
    return operator()(cell.i(), cell.j() + offset);
  case Coordinates::J_NEG:
    return operator()(cell.i(), cell.j() - offset);
  }
}

const Face &StructuredGrid2D::face(const Cell &cell,
                                   Coordinates::Direction dir) const {
  switch (dir) {
  case Coordinates::I_POS:
    return iface(cell.i() + 1, cell.j());
  case Coordinates::I_NEG:
    return iface(cell.i(), cell.j());
  case Coordinates::J_POS:
    return jface(cell.i(), cell.j() + 1);
  case Coordinates::J_NEG:
    return jface(cell.i(), cell.j());
  }
}

Scalar StructuredGrid2D::dh(const Cell &cell, Coordinates::Direction dir,
                            int offset) const {
  switch (dir) {
  case Coordinates::I_POS:
    return operator()(cell.i() + offset, cell.j()).centroid().x -
           cell.centroid().x;
  case Coordinates::I_NEG:
    return cell.centroid().x -
           operator()(cell.i() - offset, cell.j()).centroid().x;
  case Coordinates::J_POS:
    return operator()(cell.i(), cell.j() + offset).centroid().y -
           cell.centroid().y;
  case Coordinates::J_NEG:
    return cell.centroid().y -
           operator()(cell.i(), cell.j() - offset).centroid().y;
  }
}

Size StructuredGrid2D::maxInc(const Cell &cell,
                              Coordinates::Direction dir) const {
  switch (dir) {
  case Coordinates::I_POS:
    return _nCellsI - cell.i() - 1;
  case Coordinates::I_NEG:
    return cell.i();
  case Coordinates::J_POS:
    return _nCellsJ - cell.j() - 1;
  case Coordinates::J_NEG:
    return cell.j();
  }
}

std::array<Size, 2>
StructuredGrid2D::computeBlockDims(bool stripPartitioning) const {
  if (stripPartitioning)
    return std::array<Size, 2>({(Size)_comm->nProcs(), 1});

  Size i = 2, num = _comm->nProcs();

  std::deque<Size> factors(1, 1);

  while (i <= num) {
    if (num % i)
      ++i;
    else {
      factors.emplace_front(i);
      num /= i;
    }
  }

  while (factors.size() > 2) {
    std::sort(factors.begin(), factors.end());
    factors[1] *= factors[0];
    factors.pop_front();
  }

  std::array<Size, 2> result;
  result[0] = factors[0];
  result[1] = factors[1];
  return result;
}

void StructuredGrid2D::initParallel(const std::vector<int> &ownership,
                                    const std::vector<Label> &gids) {
  if (ownership.size() != _cells.size() || gids.size() != _cells.size())
    throw Exception("StructuredGrid2D", "initParallel",
                    "invalid partitioning.");

  _ownership = ownership;

  _sendBuffers.clear();
  _recvBuffers.clear();
  _sendBuffers.resize(_comm->nProcs());
  _recvBuffers.resize(_comm->nProcs());

  _globalToLocalIdMap.clear();
  _globalToLocalIdMap.reserve(gids.size());

  for (auto i = 0; i < _cells.size(); ++i) {
    _cells[i].setgid(gids[i]);
    auto insert = _globalToLocalIdMap.emplace(_cells[i].gid(), _cells[i].lid());

    if (!insert.second)
      throw Exception("StructuredGrid2D", "initParallel",
                      "duplicate global id.");

    if (_ownership[i] >= _comm->nProcs())
      throw Exception("StructuredGrid2D", "initParallel",
                      "proc " + std::to_string(_ownership[i]) +
                          " is greater than the max comm rank.");
    else if (_ownership[i] != _comm->rank()) {
      _localCells.remove(_cells[i]);
      _recvBuffers[_ownership[i]].add(_cells[i]);
    }
  }

  std::vector<std::vector<Label>> recvOrders(_comm->nProcs());
  for (auto proc = 0; proc < _comm->nProcs(); ++proc) {
    for (const Cell &c : recvBuffer(proc))
      recvOrders[proc].emplace_back(c.gid());

    _comm->isend(proc, recvOrders[proc], _comm->rank());
  }

  std::vector<std::vector<Label>> sendOrders(_comm->nProcs());
  for (auto proc = 0; proc < _comm->nProcs(); ++proc) {
    sendOrders[proc].resize(_comm->probeSize<Label>(proc, proc));
    _comm->recv(proc, sendOrders[proc], proc);
  }

  _comm->waitAll();

  for (int proc = 0; proc < _comm->nProcs(); ++proc)
    for (Label gid : sendOrders[proc]) {
      Label lid = _globalToLocalIdMap.at(gid);
      _sendBuffers[proc].add(_cells[lid]);
    }
}

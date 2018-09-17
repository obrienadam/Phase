#include <unordered_map>

#include "StructuredGrid2D.h"

StructuredGrid2D::StructuredGrid2D()
    :
      _comm(std::make_shared<Communicator>())
{
    _sendBuffers.resize(_comm->nProcs());
    _recvBuffers.resize(_comm->nProcs());
}

StructuredGrid2D::StructuredGrid2D(Size nCellsI, Size nCellsJ, Scalar lx, Scalar ly)
    :
      StructuredGrid2D()
{
    init(nCellsI, nCellsJ, std::make_pair(0., lx), std::make_pair(0., ly));
}

StructuredGrid2D::StructuredGrid2D(const Input &input)
    :
      StructuredGrid2D(
          input.caseInput().get<Size>("Grid.nCellsX"),
          input.caseInput().get<Size>("Grid.nCellsY"),
          input.caseInput().get<Scalar>("Grid.width"),
          input.caseInput().get<Scalar>("Grid.height")
          )
{

}

void StructuredGrid2D::init(Size nCellsI, Size nCellsJ, std::pair<Scalar, Scalar> xb, std::pair<Scalar, Scalar> yb)
{
    _nCellsI = nCellsI;
    _nCellsJ = nCellsJ;
    _xb = xb;
    _yb = yb;

    Scalar dx = (_xb.second - _xb.first) / _nCellsI;
    Scalar dy = (_xb.second - _xb.first) / _nCellsJ;

    _nodes.clear();

    for (auto j = 0; j < nNodesJ(); ++j)
        for (auto i = 0; i < nNodesI(); ++i)
            _nodes.push_back(Point2D(i * dx + _xb.first, j * dy + _yb.first));

    _cells.clear();

    for(auto j = 0; j < _nCellsJ; ++j)
        for(auto i = 0; i < _nCellsI; ++i)
            _cells.push_back(Cell(*this, i, j));

    _ifaces.clear();

    for(int j = 0; j < _nCellsJ; ++j)
        for(int i = 0; i < nNodesI(); ++i)
            _ifaces.push_back(Face(*this, Coordinates::I, i, j));

    _jfaces.clear();

    for(int i = 0; i < _nCellsI; ++i)
        for(int j = 0; j < nNodesJ(); ++j)
            _jfaces.push_back(Face(*this, Coordinates::J, i, j));

    _faces.clear();
    for(const Face &f: _ifaces)
        _faces.push_back(f);

    for(const Face &f: _jfaces)
        _faces.push_back(f);

    _localCells.add(_cells.begin(), _cells.end());
    _ownership.resize(_cells.size(), _comm->rank());
}

void StructuredGrid2D::init(Size nCellsI, Size nCellsJ, Scalar lx, Scalar ly, int nBufferCells)
{
    Size nBlocks[2] = {2, 2};
    Size localBlock[2] = {_comm->rank() % nBlocks[0], _comm->rank() / nBlocks[0]};

    //- map global cell ids to the owning proc
    std::vector<int> ownership(nCellsI * nCellsJ);

    auto getGlobalId = [nCellsI](Label i, Label j){ return j * nCellsI + i; };

    std::vector<Size> nLocalCellsI = _comm->allGather(nCellsI / nBlocks[0] + (nCellsI % nBlocks[0] > localBlock[0] ? 1 : 0));
    std::vector<Size> nLocalCellsJ = _comm->allGather(nCellsJ / nBlocks[1] + (nCellsJ % nBlocks[1] > localBlock[1] ? 1 : 0));

    for(auto blockJ = 0; blockJ < nBlocks[1]; ++blockJ)
        for(auto blockI = 0; blockI < nBlocks[0]; ++blockI)
        {
            int proc = blockJ * nBlocks[0] + blockI;

            std::pair<Label, Label> irange = std::make_pair(
                        std::accumulate(nLocalCellsI.begin(), nLocalCellsI.begin() + proc, 0),
                        std::accumulate(nLocalCellsI.begin(), nLocalCellsI.begin() + proc + 1, 0)
                        );

            std::pair<Label, Label> jrange = std::make_pair(
                        std::accumulate(nLocalCellsJ.begin(), nLocalCellsJ.begin() + proc, 0),
                        std::accumulate(nLocalCellsJ.begin(), nLocalCellsJ.begin() + proc + 1, 0)
                        );

            for(int j = jrange.first; j < jrange.second; ++j)
                for(int i = irange.first; i < irange.second; ++i)
                    ownership[getGlobalId(i, j)] = proc;
        }

    Size ilower = std::accumulate(nLocalCellsI.begin(), nLocalCellsI.begin() + _comm->rank(), 0);
    Size iupper = ilower + nLocalCellsI[_comm->rank()];

    Size jlower = std::accumulate(nLocalCellsJ.begin(), nLocalCellsJ.begin() + _comm->rank(), 0);
    Size jupper = jlower + nLocalCellsJ[_comm->rank()];

    ilower -= ilower > nBufferCells ? nBufferCells : ilower;
    iupper += nCellsI > iupper + nBufferCells ? nBufferCells : (nCellsI - iupper);

    jlower -= jlower > nBufferCells ? nBufferCells : jlower;
    jupper += nCellsJ > jupper + nBufferCells ? nBufferCells : (nCellsJ - jupper);

    Scalar dx = lx / nCellsI;
    Scalar dy = ly / nCellsJ;

    _nCellsI = iupper - ilower;
    _nCellsJ = jupper - jlower;
    _xb = std::make_pair(ilower * dx, iupper * dx);
    _yb = std::make_pair(jlower * dy, jupper * dy);

    init(_nCellsI, _nCellsJ, _xb, _yb);

    _ownership.clear();
    for(Cell &cell: _cells)
    {
        Label iglobal = cell.i() + ilower;
        Label jglobal = cell.j() + jlower;

        cell.setgid(getGlobalId(iglobal, jglobal));

        _ownership.push_back(ownership[cell.gid()]);

        if(_ownership.back() != _comm->rank())
        {
            _localCells.remove(cell);
            _recvBuffers[_ownership.back()].add(cell);
        }
    }

    std::unordered_map<Label, Label> globalToLocalIdMap;
    for(const Cell &cell: _cells)
        globalToLocalIdMap[cell.gid()] = cell.lid();

    //- Send the correct recv order to neighbouring procs
    std::vector<std::vector<Label>> recvOrders(_comm->nProcs());
    for(int proc = 0; proc < _comm->nProcs(); ++proc)
    {
        recvOrders[proc].resize(_recvBuffers[proc].size());

        std::transform(_recvBuffers[proc].begin(),
                       _recvBuffers[proc].end(),
                       recvOrders[proc].begin(),
                       [](const Cell &c) { return c.gid(); });

        _comm->isend(proc, recvOrders[proc], _comm->rank());
    }

    //- Initialize the send buffers
    std::vector<Label> sendOrder;
    for(int proc = 0; proc < _comm->nProcs(); ++proc)
    {
        //- Get the global ids
        sendOrder.resize(_comm->probeSize<Label>(proc, proc));
        _comm->recv(proc, sendOrder, proc);

        //- map global ids to local ids
        std::transform(sendOrder.begin(),
                       sendOrder.end(),
                       sendOrder.begin(),
                       [&globalToLocalIdMap](Label gid) { return globalToLocalIdMap.at(gid); });

        for(Label lid: sendOrder)
            _sendBuffers[proc].add(_cells[lid]);
    }
}

const Cell &StructuredGrid2D::cell(const Cell &cell, Coordinates::Direction dir, int offset) const
{
    switch(dir)
    {
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

const Face &StructuredGrid2D::face(const Cell &cell, Coordinates::Direction dir) const
{
    switch(dir)
    {
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

Scalar StructuredGrid2D::dh(const Cell &cell, Coordinates::Direction dir, int offset) const
{
    switch(dir)
    {
    case Coordinates::I_POS:
        return operator()(cell.i() + offset, cell.j()).centroid().x - cell.centroid().x;
    case Coordinates::I_NEG:
        return cell.centroid().x - operator()(cell.i() - offset, cell.j()).centroid().x;
    case Coordinates::J_POS:
        return operator()(cell.i(), cell.j() + offset).centroid().y - cell.centroid().y;
    case Coordinates::J_NEG:
        return cell.centroid().y - operator()(cell.i(), cell.j() - offset).centroid().y;
    }
}

Size StructuredGrid2D::maxInc(const Cell &cell, Coordinates::Direction dir) const
{
    switch(dir)
    {
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

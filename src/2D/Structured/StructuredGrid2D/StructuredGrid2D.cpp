#include "StructuredGrid2D.h"

StructuredGrid2D::StructuredGrid2D(Size nCellsI, Size nCellsJ, Scalar lx, Scalar ly)
    :
      StructuredGrid2D()
{
    init(nCellsI, nCellsJ, lx, ly);
}

StructuredGrid2D::StructuredGrid2D(const std::vector<Scalar> &xcoords, const std::vector<Scalar> &ycoords)
    :
      StructuredGrid2D()
{
    init(xcoords, ycoords);
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

void StructuredGrid2D::init(Size nCellsI, Size nCellsJ, Scalar lx, Scalar ly)
{
    _nCellsI = nCellsI;
    _nCellsJ = nCellsJ;
    _lx = lx;
    _ly = ly;

    Scalar dx = _lx / _nCellsI;
    Scalar dy = _ly / _nCellsJ;

    _nodes.clear();

    for (auto j = 0; j < nNodesJ(); ++j)
        for (auto i = 0; i < nNodesI(); ++i)
            _nodes.push_back(Point2D(i * dx, j * dy));

    _cells.clear();

    for(auto j = 0; j < _nCellsJ; ++j)
        for(auto i = 0; i < _nCellsI; ++i)
            _cells.push_back(Cell(*this, i, j));

    _ifaces.clear();

    for(int j = 0; j < _nCellsJ; ++j)
        for(int i = 0; i < nNodesI(); ++i)
            _ifaces.push_back(Face(*this, I, i, j));

    _jfaces.clear();

    for(int i = 0; i < _nCellsI; ++i)
        for(int j = 0; j < nNodesJ(); ++j)
            _jfaces.push_back(Face(*this, J, i, j));

    _faces.clear();
    for(const Face &f: _ifaces)
        _faces.push_back(f);

    for(const Face &f: _jfaces)
        _faces.push_back(f);

    _localCells.add(_cells.begin(), _cells.end());
    _ownership.resize(_cells.size(), 0);
}

void StructuredGrid2D::init(std::vector<Scalar> xcoords, std::vector<Scalar> ycoords)
{
    std::sort(xcoords.begin(), xcoords.end());
    std::sort(ycoords.begin(), ycoords.end());

    _nCellsI = xcoords.size() - 1;
    _nCellsJ = ycoords.size() - 1;
    _lx = xcoords.back();
    _ly = ycoords.back();

    _nodes.clear();

    for(auto j = 0; j < nNodesJ(); ++j)
        for(auto i = 0; i < nNodesI(); ++i)
            _nodes.push_back(Point2D(xcoords[i], ycoords[j]));

    _cells.clear();

    for(auto j = 0; j < _nCellsJ; ++j)
        for(auto i = 0; i < _nCellsI; ++i)
            _cells.push_back(Cell(*this, i, j));

    _ifaces.clear();

    for(int j = 0; j < _nCellsJ; ++j)
        for(int i = 0; i < nNodesI(); ++i)
            _ifaces.push_back(Face(*this, I, i, j));

    _jfaces.clear();

    for(int i = 0; i < _nCellsI; ++i)
        for(int j = 0; j < nNodesJ(); ++j)
            _jfaces.push_back(Face(*this, J, i, j));

    _faces.clear();
    for(const Face &f: _ifaces)
        _faces.push_back(f);

    for(const Face &f: _jfaces)
        _faces.push_back(f);

    _localCells.add(_cells.begin(), _cells.end());
    _ownership.resize(_cells.size(), 0);
}

const Cell &StructuredGrid2D::cell(const Cell &cell, Orientation dir, int offset) const
{
    switch(dir)
    {
    case I_POS:
        return operator()(cell.i() + offset, cell.j());
    case I_NEG:
        return operator()(cell.i() - offset, cell.j());
    case J_POS:
        return operator()(cell.i(), cell.j() + offset);
    case J_NEG:
        return operator()(cell.i(), cell.j() - offset);
    }
}

const Face &StructuredGrid2D::face(const Cell &cell, Orientation dir) const
{
    switch(dir)
    {
    case I_POS:
        return iface(cell.i() + 1, cell.j());
    case I_NEG:
        return iface(cell.i(), cell.j());
    case J_POS:
        return jface(cell.i(), cell.j() + 1);
    case J_NEG:
        return jface(cell.i(), cell.j());
    }
}

Scalar StructuredGrid2D::dh(const Cell &cell, Orientation dir, int offset) const
{
    switch(dir)
    {
    case I_POS:
        return operator()(cell.i() + offset, cell.j()).centroid().x - cell.centroid().x;
    case I_NEG:
        return cell.centroid().x - operator()(cell.i() - offset, cell.j()).centroid().x;
    case J_POS:
        return operator()(cell.i(), cell.j() + offset).centroid().y - cell.centroid().y;
    case J_NEG:
        return cell.centroid().y - operator()(cell.i(), cell.j() - offset).centroid().y;
    }
}

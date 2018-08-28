#include "StructuredGrid2D.h"

StructuredGrid2D::StructuredGrid2D(Size nCellsI, Size nCellsJ, Scalar lx, Scalar ly)
{
    init(nCellsI, nCellsJ, lx, ly);
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
            _ifaces.push_back(Face(*this, Face::I, i, j));

    _jfaces.clear();

    for(int i = 0; i < _nCellsI; ++i)
        for(int j = 0; j < nNodesJ(); ++j)
            _jfaces.push_back(Face(*this, Face::J, i, j));

    _faces.clear();
    for(const Face &f: _ifaces)
        _faces.push_back(f);

    for(const Face &f: _jfaces)
        _faces.push_back(f);
}

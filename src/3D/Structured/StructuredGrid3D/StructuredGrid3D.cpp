#include "StructuredGrid3D.h"

StructuredGrid3D::StructuredGrid3D(Size nCellsI, Size nCellsJ, Size nCellsK, Scalar lx, Scalar ly, Scalar lz)
    :
      _nCellsI(nCellsI),
      _nCellsJ(nCellsJ),
      _nCellsK(nCellsK),
      _nCells(nCellsI * nCellsJ * nCellsK),
      _lx(lx),
      _ly(ly),
      _lz(lz)
{
    Scalar dx = _lx / _nCellsI;
    Scalar dy = _ly / _nCellsJ;
    Scalar dz = _lz / _nCellsK;

    for(int k = 0; k < _nCellsK + 1; ++k)
        for(int j = 0; j < _nCellsJ + 1; ++j)
            for(int i = 0; i < _nCellsI + 1; ++i)
                _nodes.push_back(Point3D(i * dx, j * dy, k * dz));

    for(int k = 0; k < _nCellsK; ++k)
        for(int j = 0; j < _nCellsJ; ++j)
            for(int i = 0; i < _nCellsI; ++i)
            {
                _cells.push_back(Cell(*this, i, j, k));

                //- Initialize faces
                if(i == 0)
                    _ifaces.push_back(Face(*this, Face::I_NEG, i, j, k));

                _ifaces.push_back(Face(*this, Face::I_POS, i, j, k));

                if(j == 0)
                    _jfaces.push_back(Face(*this, Face::J_NEG, i, j, k));

                _jfaces.push_back(Face(*this, Face::J_POS, i, j, k));

                if(k == 0)
                    _kfaces.push_back(Face(*this, Face::K_NEG, i, j, k));

                _kfaces.push_back(Face(*this, Face::K_POS, i, j, k));
            }

    _comm = std::make_shared<Communicator>();
}

StructuredGrid3D::StructuredGrid3D(const Input &input)
    :
      StructuredGrid3D(
          input.caseInput().get<Size>("Grid.nCellsX"),
          input.caseInput().get<Size>("Grid.nCellsY"),
          input.caseInput().get<Size>("Grid.nCellsZ"),
          input.caseInput().get<Scalar>("Grid.width"),
          input.caseInput().get<Scalar>("Grid.height"),
          input.caseInput().get<Scalar>("Grid.depth")
          )
{

}

const Face &StructuredGrid3D::operator()(Label i, Label j, Label k, Face::Direction f) const
{
    switch (f)
    {
    case Face::I_POS:
        return _ifaces[k * (_nCellsI + 1) * _nCellsJ + j * (_nCellsI + 1) + i + 1];
    case Face::I_NEG:
        return _ifaces[k * (_nCellsI + 1) * _nCellsJ + j * (_nCellsI + 1) + i];
    case Face::J_POS:
        return _jfaces[k * _nCellsI * (_nCellsJ + 1) + (j + 1) * _nCellsI + i];
    case Face::J_NEG:
        return _jfaces[k * _nCellsI * (_nCellsJ + 1) + j * _nCellsI + i];
    case Face::K_POS:
        return _kfaces[(k + 1) * _nCellsI * _nCellsJ + j * _nCellsI + i];
    case Face::K_NEG:
        return _kfaces[k * _nCellsI * _nCellsJ + j * _nCellsI + i];
    }
}

const Cell &StructuredGrid3D::operator()(Label i, Label j, Label k) const
{
    return _cells[_nCellsI * _nCellsJ * k + _nCellsI * j + i];
}

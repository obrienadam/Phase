#include "FaceStencil.h"
#include "Cell.h"
#include "StructuredGrid3D.h"

FaceStencil::FaceStencil(const Cell &cell, Face::Direction dir, int forwardBias, int backwardBias)
    :
      _cell(cell),
      _face(cell.grid()(cell, dir)),
      _dir(dir)
{
    //- Shift the stencil if necessary
    if(forwardBias > maxForwardShift())
    {
        int shift = forwardBias - maxForwardShift();
        forwardBias -= shift;
        backwardBias += shift + 1;
    }
    else if(backwardBias > maxBackwardShift() + 1)
    {
        int shift = backwardBias - maxBackwardShift();
        forwardBias += shift + 1;
        backwardBias -= shift;
    }

    _cells.reserve(forwardBias + backwardBias);
    for(int i = -backwardBias + 1; i <= forwardBias; ++i)
        _cells.push_back(std::cref(this->operator()(i)));

    std::vector<Scalar> offsets;
    offsets.reserve(_cells.size() + 1);

    for(const Cell& stCell: _cells)
    {
        switch(_dir)
        {
        case Face::I_POS:
            offsets.push_back(stCell.centroid().x - _face.centroid().x);
            break;
        case Face::I_NEG:
            offsets.push_back(_face.centroid().x - stCell.centroid().x);
            break;
        case Face::J_POS:
            offsets.push_back(stCell.centroid().y - _face.centroid().y);
            break;
        case Face::J_NEG:
            offsets.push_back(_face.centroid().y - stCell.centroid().y);
            break;
        case Face::K_POS:
            offsets.push_back(stCell.centroid().z - _face.centroid().z);
            break;
        case Face::K_NEG:
            offsets.push_back(_face.centroid().z - stCell.centroid().z);
            break;
        }
    }

    //- The face location
    offsets.push_back(0.);

    _taylorCoeffs.computeTaylorCoeffs(offsets, 1);
}

FaceStencil::FaceStencil(const Cell &cell, Face::Direction dir, int order)
    :
      FaceStencil(cell, dir, order / 2, order / 2)
{

}

Size FaceStencil::maxForwardShift()
{
    switch (_dir)
    {
    case Face::I_POS: return _cell.grid().iMax() - _cell.i();
    case Face::I_NEG: return _cell.i();
    case Face::J_POS: return _cell.grid().jMax() - _cell.j();
    case Face::J_NEG: return _cell.j();
    case Face::K_POS: return _cell.grid().kMax() - _cell.k();
    case Face::K_NEG: return _cell.k();
    }
}

Size FaceStencil::maxBackwardShift()
{
    switch (_dir)
    {
    case Face::I_POS: return _cell.i();
    case Face::I_NEG: return _cell.grid().iMax() - _cell.i();
    case Face::J_POS: return _cell.j();
    case Face::J_NEG: return _cell.grid().jMax() - _cell.j();
    case Face::K_POS: return _cell.k();
    case Face::K_NEG: return _cell.grid().kMax() - _cell.k();
    }
}

const Cell &FaceStencil::operator()(Size shift) const
{
    Label i = _cell.i(), j = _cell.j(), k = _cell.k();

    switch (_dir)
    {
    case Face::I_POS: return _cell.grid()(i + shift, j, k);
    case Face::I_NEG: return _cell.grid()(i - shift, j, k);
    case Face::J_POS: return _cell.grid()(i, j + shift, k);
    case Face::J_NEG: return _cell.grid()(i, j - shift, k);
    case Face::K_POS: return _cell.grid()(i, j, k + shift);
    case Face::K_NEG: return _cell.grid()(i, j, k - shift);
    }
}

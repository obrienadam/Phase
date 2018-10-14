#include "DirectForcingImmersedBoundaryLeastSquaresQuadraticStencil.h"

std::vector<std::vector<const ImmersedBoundaryObject*>> DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil::_ibObjSets(2);

Matrix DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil::_A;

Matrix DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil::_b;

DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil::LeastSquaresQuadraticStencil(const Cell &cell,
                                                                                          const DirectForcingImmersedBoundary &ib)
{
    _ibObjSets[0].clear();
    for(const CellLink &nb: cell.neighbours())
    {
        auto ibObj = ib.ibObj(nb.cell());

        if(ibObj && std::find(_ibObjSets[0].begin(), _ibObjSets[0].end(), ibObj.get()) == _ibObjSets[0].end())
        {
            _compatPts.push_back(CompatPoint(cell, *ibObj));
            _ibObjSets[0].emplace_back(ibObj.get());
        }
        else
            _cells.push_back(&nb.cell());
    }

    for(const CellLink &nb: cell.diagonals())
    {
        auto ibObj = ib.ibObj(nb.cell());

        if(!ibObj)
            _cells.push_back(&nb.cell());
    }

    for(const BoundaryLink &bd: cell.boundaries())
    {
        auto ibObj = ib.ibObj(bd.face().centroid());
        if(!ibObj)
            _faces.push_back(&bd.face());
    }

    for(const Cell *stCell: _cells)
    {
        _ibObjSets[1].clear();

        for(const CellLink &nb: stCell->neighbours())
        {
            auto ibObj = ib.ibObj(nb.cell());

            if(ibObj && std::find(_ibObjSets[0].begin(), _ibObjSets[0].end(), ibObj.get()) != _ibObjSets[0].end()
                    && std::find(_ibObjSets[1].begin(), _ibObjSets[1].end(), ibObj.get()) == _ibObjSets[1].end())
            {
                _compatPts.push_back(CompatPoint(*stCell, *ibObj));
                _ibObjSets[1].emplace_back(ibObj.get());
            }
        }
    }

    if(nReconstructionPoints() < 2)
        throw Exception("DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil",
                        "LeastSquaresQuadraticStencil",
                        "Could not locate enough reconstruction points. Points found = "
                        + std::to_string(nReconstructionPoints()) + "."
                        + " Num cells = " + std::to_string(_cells.size())
                        + ", Num compat pts = " + std::to_string(_compatPts.size()) + ".");
}

Matrix DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil::interpolationCoeffs(const Point2D &x) const
{
    if(nReconstructionPoints() >= 6)
        return quadraticInterpolationCoeffs(x);
    else if(nReconstructionPoints() >= 3)
        return linearInterpolationCoeffs(x);
    else if(_compatPts.size()  == 2)
        return subgridInterpolationCoeffs(x);
}

Matrix DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil::linearInterpolationCoeffs(const Point2D &x) const
{
    _A.resize(nReconstructionPoints(), 3);

    int i = 0;
    for(const Cell *cell: _cells)
    {
        const Point2D &x = cell->centroid();
        _A.setRow(i++, {x.x, x.y, 1.});
    }

    for(const CompatPoint &cpt: _compatPts)
    {
        const Point2D &x = cpt.pt();
        _A.setRow(i++, {x.x, x.y, 1.});
    }

    for(const Face *face: _faces)
    {
        const Point2D &x = face->centroid();
        _A.setRow(i++, {x.x, x.y, 1.});
    }

    _b.resize(1, 3);
    _b.setRow(0, {x.x, x.y, 1.});
    _A.pinvert();

    return _b * _A;
}

Matrix DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil::quadraticInterpolationCoeffs(const Point2D &x) const
{
    _A.resize(nReconstructionPoints(), 6);

    int i = 0;
    for(const Cell *cell: _cells)
    {
        const Point2D &x = cell->centroid();
        _A.setRow(i++, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});
    }

    for(const CompatPoint &cpt: _compatPts)
    {
        const Point2D &x = cpt.pt();
        _A.setRow(i++, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});
    }

    for(const Face *face: _faces)
    {
        const Point2D &x = face->centroid();
        _A.setRow(i++, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});
    }

    _b.resize(1, 6);
    _b.setRow(0, {x.x * x.x, x.y * x.y, x.x * x.y, x.x, x.y, 1.});
    _A.pinvert();

    return _b * _A;
}

Matrix DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil::subgridInterpolationCoeffs(const Point2D &x) const
{
    Point2D pt1 = _compatPts[0].pt();
    Point2D pt2 = _compatPts[1].pt();

    Scalar l1 = (pt1 - x).mag();
    Scalar l2 = (pt2 - x).mag();

    Scalar g = l2 / (l1 + l2);

    return Matrix(1, 2, {g, 1. - g});
}

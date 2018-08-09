#ifndef PHASE_DIRECT_FORCING_IMMERSED_BOUNDARY_LEAST_SQUARES_QUADRATIC_STENCIL_H
#define PHASE_DIRECT_FORCING_IMMERSED_BOUNDARY_LEAST_SQUARES_QUADRATIC_STENCIL_H

#include "DirectForcingImmersedBoundary.h"

class DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil
{
public:

    class CompatPoint
    {
    public:

        CompatPoint(const Cell &cell, const ImmersedBoundaryObject &ibObj)
            : _cell(&cell), _ibObj(&ibObj), _pt(ibObj.nearestIntersect(cell.centroid()))
        {}

        Vector2D velocity() const
        { return _ibObj->velocity(_pt); }

        Vector2D acceleration() const
        { return _ibObj->acceleration(_pt); }

        const Cell &cell() const
        { return *_cell; }

        const ImmersedBoundaryObject &ibObj() const
        { return *_ibObj; }

        const Point2D &pt() const
        { return _pt; }

    private:

        const Cell *_cell;

        const ImmersedBoundaryObject *_ibObj;

        Point2D _pt;
    };

    LeastSquaresQuadraticStencil(const Cell &cell,
                                 const DirectForcingImmersedBoundary &ib);

    Size nReconstructionPoints() const
    { return _cells.size() + _compatPts.size(); }

    const std::vector<const Cell*> &cells() const
    { return _cells; }

    const std::vector<CompatPoint> &compatPts() const
    { return _compatPts; }

protected:

    std::vector<const Cell*> _cells;

    std::vector<CompatPoint> _compatPts;
};

#endif

#ifndef PHASE_DIRECT_FORCING_IMMERSED_BOUNDARY_LEAST_SQUARES_QUADRATIC_STENCIL_H
#define PHASE_DIRECT_FORCING_IMMERSED_BOUNDARY_LEAST_SQUARES_QUADRATIC_STENCIL_H

#include "System/StaticVector.h"

#include "DirectForcingImmersedBoundary.h"

class DirectForcingImmersedBoundary::LeastSquaresQuadraticStencil
{
public:

    class CompatPoint
    {
    public:

        CompatPoint()
        {}

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

        Vector2D ns() const
        { return _ibObj->nearestEdgeUnitNormal(_pt); }

    private:

        const Cell *_cell;

        const ImmersedBoundaryObject *_ibObj;

        Point2D _pt;
    };

    LeastSquaresQuadraticStencil(const Cell &cell,
                                 const DirectForcingImmersedBoundary &ib);

    Size nReconstructionPoints() const
    { return _cells.size() + _faces.size() + _compatPts.size(); }

    const StaticVector<const Cell*, 8> &cells() const
    { return _cells; }

    const StaticVector<const Face*, 8> &faces() const
    { return _faces; }

    const StaticVector<CompatPoint, 8> &compatPts() const
    { return _compatPts; }

    //
    Matrix interpolationCoeffs(const Point2D &x) const;

    Matrix continuityConstrainedInterpolationCoeffs(const Point2D &pt) const;

    Matrix polarQuadraticContinuityConstrainedInterpolationCoeffs(const Point2D &x) const;

protected:

    static std::vector<std::vector<const ImmersedBoundaryObject*>> _ibObjSets;

    static Matrix _A, _b;

    Matrix linearInterpolationCoeffs(const Point2D &x) const;

    Matrix quadraticInterpolationCoeffs(const Point2D &x) const;

    Matrix subgridInterpolationCoeffs(const Point2D &x) const;

    Matrix quadraticContinuityConstrainedInterpolationCoeffs(const Point2D &x) const;

    StaticVector<const Cell*, 8> _cells;

    StaticVector<const Face*, 8> _faces;

    StaticVector<CompatPoint, 8> _compatPts;
};

#endif

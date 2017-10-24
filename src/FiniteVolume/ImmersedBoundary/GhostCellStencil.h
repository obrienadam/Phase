#ifndef GHOST_CELL_STENCIL_H
#define GHOST_CELL_STENCIL_H

#include "Cell.h"
#include "Interpolation.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"
#include "ImmersedBoundaryStencil.h"
#include "StaticMatrix.h"

class GhostCellImmersedBoundaryObject;

class GhostCellStencil: public ImmersedBoundaryStencil
{
public:

    GhostCellStencil(const Cell& cell): ImmersedBoundaryStencil(cell) {}

    GhostCellStencil(const Cell &cell, const GhostCellImmersedBoundaryObject &ibObj, const FiniteVolumeGrid2D &grid);

    GhostCellStencil(const Cell& cell,
                     const Point2D& bp,
                     const Vector2D& cl,
                     const FiniteVolumeGrid2D &grid);

    const Point2D &boundaryPoint() const
    { return bp_; }

    const Point2D &imagePoint() const
    { return ip_; }

    const std::vector<Ref<const Cell>>& cells() const
    { return cells_; }

    Vector2D unitNormal() const
    { return (bp_ - ip_).unitVec(); }

    Scalar length() const
    { return (ip_ - cell().centroid()).mag(); }

    bool diagonalStencil() const;

    Scalar ipValue(const ScalarFiniteVolumeField& field) const;

    Vector2D ipValue(const VectorFiniteVolumeField& field) const;

    Scalar bpValue(const ScalarFiniteVolumeField& field) const;

    Vector2D bpValue(const VectorFiniteVolumeField& field) const;

    Vector2D ipGrad(const ScalarFiniteVolumeField& field) const;

    Vector2D bpGrad(const ScalarFiniteVolumeField& field) const;

protected:

    StaticMatrix<4, 4> A_;
    Point2D ip_, bp_;
};

#endif

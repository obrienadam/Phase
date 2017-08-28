#ifndef GHOST_CELL_STENCIL_H
#define GHOST_CELL_STENCIL_H

#include "Cell.h"
#include "Interpolation.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"
#include "ImmersedBoundaryStencil.h"

class GhostCellImmersedBoundaryObject;

class GhostCellStencil: public ImmersedBoundaryStencil
{
public:

    GhostCellStencil(const Cell &cell, const GhostCellImmersedBoundaryObject &ibObj, const FiniteVolumeGrid2D &grid);

    GhostCellStencil(const Cell &cell,
                     const GhostCellImmersedBoundaryObject &ibObj,
                     const FiniteVolumeGrid2D &grid,
                     Scalar theta,
                     const Vector2D& m);

    const Point2D &boundaryPoint() const
    { return bp_; }

    const std::vector<Ref<const Cell>>& cells() const
    { return cells_; }

    const std::vector<Scalar>& dirichletCoeffs() const
    { return dirichletCoeffs_; }

    const std::vector<Scalar>& neumannCoeffs() const
    { return neumannCoeffs_; }

    Vector2D unitNormal() const
    { return (bp_ - ip_).unitVec(); }

    Scalar length() const
    { return (ip_ - cell_.centroid()).mag(); }

    Scalar ipValue(const ScalarFiniteVolumeField& field) const;

    Vector2D ipValue(const VectorFiniteVolumeField& field) const;

    Vector2D ipGrad(const ScalarFiniteVolumeField& field) const;

protected:

    Matrix A_;

    Point2D ip_, bp_;
    std::vector<Ref<const Cell>> cells_;
    std::vector<Scalar> dirichletCoeffs_, neumannCoeffs_;
};

#endif

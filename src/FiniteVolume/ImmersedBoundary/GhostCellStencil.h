#ifndef GHOST_CELL_STENCIL_H
#define GHOST_CELL_STENCIL_H

#include "ImmersedBoundaryStencil.h"
#include "Cell.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"
#include "StaticMatrix.h"

class GhostCellStencil : public ImmersedBoundaryStencil
{
public:

    GhostCellStencil(const Cell &cell) : ImmersedBoundaryStencil(cell)
    {}

    GhostCellStencil(const Cell &cell,
                     const ImmersedBoundaryObject &ibObj,
                     bool throwExceptionOnError = true);

    GhostCellStencil(const Cell &cell,
                     const Point2D &bp,
                     const Vector2D &cl,
                     const FiniteVolumeGrid2D &grid,
                     bool throwExceptionOnError = true);


    const Point2D &boundaryPoint() const
    { return bp_; }

    const Point2D &imagePoint() const
    { return ip_; }

    const Vector2D &wallNormal() const
    {
        return nw_;
    }

    const std::vector<Ref<const Cell>> &cells() const
    { return cells_; }

    Scalar length() const
    { return (ip_ - cell().centroid()).mag(); }

    Scalar ipValue(const ScalarFiniteVolumeField &field) const;

    Vector2D ipValue(const VectorFiniteVolumeField &field) const;

    Scalar bpValue(const ScalarFiniteVolumeField &field) const;

    Vector2D bpValue(const VectorFiniteVolumeField &field) const;

    Vector2D ipGrad(const ScalarFiniteVolumeField &field) const;

    Vector2D bpGrad(const ScalarFiniteVolumeField &field) const;

    Tensor2D bpGrad(const VectorFiniteVolumeField &field) const;

protected:

    StaticMatrix<4, 4> A_;
    Point2D ip_, bp_, nw_;
};

#endif

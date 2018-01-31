#ifndef GHOST_CELL_STENCIL_H
#define GHOST_CELL_STENCIL_H

#include "ImmersedBoundaryStencil.h"
#include "Cell.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"
#include "StaticMatrix.h"

class GhostCellStencil
{
public:

    GhostCellStencil(const Cell &cell,
                     const ImmersedBoundaryObject &ibObj);

    GhostCellStencil(const Cell &cell,
                     const Point2D &bp,
                     const Vector2D &cl,
                     const FiniteVolumeGrid2D &grid);

    const Cell &cell() const
    { return cell_; }

    const std::vector<Ref<const Cell>> &dirichletCells() const
    { return dirichletCells_; }

    const std::vector<Ref<const Cell>> &neumannCells() const
    { return neumannCells_; }

    const std::vector<Scalar> &dirichletCoeffs() const
    { return dirichletCoeffs_; }

    const std::vector<Scalar> &neumannCoeffs() const
    { return neumannCoeffs_; }

    bool ghostCellInDirichletStencil() const
    { return ghostCellInDirichletStencil_; }

    bool ghostCellInNeumannStencil() const
    { return ghostCellInNeumannStencil_; }

    const Point2D &boundaryPoint() const
    { return bp_; }

    const Point2D &imagePoint() const
    { return ip_; }

    const Vector2D &wallNormal() const
    { return nw_; }

    Scalar length() const
    { return (ip_ - cell_.get().centroid()).mag(); }

    Scalar ipValue(const ScalarFiniteVolumeField &field) const;

    Vector2D ipValue(const VectorFiniteVolumeField &field) const;

    Scalar bpValue(const ScalarFiniteVolumeField &field) const;

    Vector2D bpValue(const VectorFiniteVolumeField &field) const;

    Vector2D ipGrad(const ScalarFiniteVolumeField &field) const;

    Vector2D bpGrad(const ScalarFiniteVolumeField &field) const;

    Tensor2D bpGrad(const VectorFiniteVolumeField &field) const;

protected:

    void initDirichletCoeffs();

    void initNeumannCoeffs();

    StaticMatrix<4, 4> Ad_, An_;

    Ref<const Cell> cell_;

    std::vector<Ref<const Cell>> dirichletCells_, neumannCells_;

    std::vector<Scalar> dirichletCoeffs_, neumannCoeffs_;

    bool ghostCellInNeumannStencil_, ghostCellInDirichletStencil_;

    Point2D ip_, bp_, nw_;
};

#endif

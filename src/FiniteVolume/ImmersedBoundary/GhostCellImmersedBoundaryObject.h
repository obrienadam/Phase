#ifndef GHOST_CELL_IMMERSED_BOUNDARY_OBJECT_H
#define GHOST_CELL_IMMERSED_BOUNDARY_OBJECT_H

#include "ImmersedBoundaryObject.h"
#include "GhostCellStencil.h"
#include "Multiphase/Celeste.h"

class GhostCellImmersedBoundaryObject : public ImmersedBoundaryObject
{
public:

    GhostCellImmersedBoundaryObject(const std::string &name,
                                    Label id,
                                    const ImmersedBoundary &ib,
                                    const std::shared_ptr<FiniteVolumeGrid2D>& grid);

    Type type() const
    { return GHOST_CELL; }

    void updateCells();

    void updateContactLineStencils(Scalar theta);

    Equation<Scalar> bcs(ScalarFiniteVolumeField &field) const;

    Equation<Vector2D> bcs(VectorFiniteVolumeField &field) const;

    Equation<Vector2D> velocityBcs(VectorFiniteVolumeField &u) const;

    Equation<Scalar> pressureBcs(Scalar rho, ScalarFiniteVolumeField &p) const;

    Equation<Scalar> contactLineBcs(ScalarFiniteVolumeField &gamma, Scalar theta) const;

    void computeForce(Scalar rho,
                      Scalar mu,
                      const VectorFiniteVolumeField &u,
                      const ScalarFiniteVolumeField &p,
                      const Vector2D &g = Vector2D(0., 0.));

    void computeForce(const ScalarFiniteVolumeField &rho,
                      const ScalarFiniteVolumeField &mu,
                      const VectorFiniteVolumeField &u,
                      const ScalarFiniteVolumeField &p,
                      const Vector2D &g = Vector2D(0., 0.));


    const std::vector<GhostCellStencil> &stencils() const
    { return stencils_; }

    const std::vector<GhostCellStencil> &contactLineStencils() const
    { return contactLineStencils_; }

protected:

    void constructStencils();

    std::vector<GhostCellStencil> stencils_;
    std::vector<GhostCellStencil> contactLineStencils_;
};

#endif

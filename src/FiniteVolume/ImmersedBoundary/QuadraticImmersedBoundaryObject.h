#ifndef QUADRATIC_IMMERSED_BOUNDARY_OBJECT_H
#define QUADRATIC_IMMERSED_BOUNDARY_OBJECT_H

#include "ImmersedBoundaryObject.h"
#include "QuadraticIbmStencil.h"

class QuadraticImmersedBoundaryObject : public ImmersedBoundaryObject
{
public:

    //- Constructors, one for circles, another for polygons
    QuadraticImmersedBoundaryObject(const std::string &name,
                                    Label id,
                                    const ImmersedBoundary &ib,
                                    const std::shared_ptr<FiniteVolumeGrid2D> &grid);

    Type type() const
    { return QUADRATIC; }

    //- Clear

    void clear();

    //- Update
    void updateCells();

    //- Boundary conditions
    Equation<Scalar> bcs(ScalarFiniteVolumeField &field) const;

    Equation<Vector2D> bcs(VectorFiniteVolumeField &field) const;

    Equation<Vector2D> velocityBcs(VectorFiniteVolumeField &u) const;

    //- Forcing

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

    void computeForce(const ScalarFiniteVolumeField &rho,
                      const ScalarFiniteVolumeField &mu,
                      const VectorFiniteVolumeField &u,
                      const ScalarFiniteVolumeField &p,
                      const ScalarFiniteVolumeField &gamma,
                      const SurfaceTensionForce &ft,
                      const Vector2D &g = Vector2D(0., 0.));

    const CellGroup &forcingCells() const
    { return forcingCells_; }

    const std::vector<QuadraticIbmStencil> &stencils() const
    { return stencils_; }

protected:

    CellGroup forcingCells_;
    std::vector<QuadraticIbmStencil> stencils_;

};

#endif

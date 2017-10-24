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
                           FiniteVolumeGrid2D &grid);

    Type type() const
    { return QUADRATIC; }

    //- Update
    void updateCells();

    //- Boundary conditions
    Equation<Scalar> bcs(ScalarFiniteVolumeField &field) const;

    Equation<Vector2D> bcs(VectorFiniteVolumeField &field) const;

    Equation<Vector2D> velocityBcs(VectorFiniteVolumeField &u) const;

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

    const CellGroup& forcingCells() const
    { return forcingCells_; }

protected:

    CellGroup forcingCells_;

};

#endif

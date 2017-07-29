#ifndef STEP_IMMERSED_BOUNDARY_OBJECT
#define STEP_IMMERSED_BOUNDARY_OBJECT

#include "ImmersedBoundaryObject.h"

class StepImmersedBoundaryObject: public ImmersedBoundaryObject
{
public:

    using ImmersedBoundaryObject::ImmersedBoundaryObject;

    virtual void update(Scalar timeStep);

    virtual void updateCells();

    virtual Equation<Scalar> bcs(ScalarFiniteVolumeField& field) const;

    virtual Equation<Vector2D> bcs(VectorFiniteVolumeField& field) const;

    void computeForce(Scalar rho, Scalar mu, const VectorFiniteVolumeField& u, const ScalarFiniteVolumeField& p);

private:
};

#endif
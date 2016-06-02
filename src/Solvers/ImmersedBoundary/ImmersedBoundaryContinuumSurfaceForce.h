#ifndef IMMERSED_BOUNDARY_CONTINUUM_SURFACE_FORCE_H
#define IMMERSED_BOUNDARY_CONTINUUM_SURFACE_FORCE_H

#include "ContinuumSurfaceForce.h"

class ImmersedBoundaryObject;

class ImmersedBoundaryContinuumSurfaceForce : public ContinuumSurfaceForce
{
public:
    ImmersedBoundaryContinuumSurfaceForce(const Input& input, const ScalarFiniteVolumeField &gamma, const ScalarFiniteVolumeField& rho);

    virtual VectorFiniteVolumeField compute(const ImmersedBoundaryObject& ibObj);

protected:

    virtual void computeInterfaceNormals(const ImmersedBoundaryObject &ibObj);

};

#endif

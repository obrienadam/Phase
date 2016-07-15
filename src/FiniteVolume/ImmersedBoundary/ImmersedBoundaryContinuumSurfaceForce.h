#ifndef IMMERSED_BOUNDARY_CONTINUUM_SURFACE_FORCE_H
#define IMMERSED_BOUNDARY_CONTINUUM_SURFACE_FORCE_H

#include "Celeste.h"

class ImmersedBoundaryObject;

class ImmersedBoundaryContinuumSurfaceForce : public ContinuumSurfaceForce
{
public:
    ImmersedBoundaryContinuumSurfaceForce(const Input& input,
                                          Solver &solver);

    virtual VectorFiniteVolumeField compute(const std::vector<ImmersedBoundaryObject> &ibObjs);

protected:

    virtual void computeInterfaceNormals(const std::vector<ImmersedBoundaryObject> &ibObjs);

};

#endif

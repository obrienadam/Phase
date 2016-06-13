#ifndef IMMERSED_BOUNDARY_CONTINUUM_SURFACE_FORCE_H
#define IMMERSED_BOUNDARY_CONTINUUM_SURFACE_FORCE_H

#include "ContinuumSurfaceForce.h"

class ImmersedBoundaryObject;

class ImmersedBoundaryContinuumSurfaceForce : public ContinuumSurfaceForce
{
public:
    ImmersedBoundaryContinuumSurfaceForce(const Input& input, const ScalarFiniteVolumeField &gamma, const VectorFiniteVolumeField& u, std::map<std::string, ScalarFiniteVolumeField>& fields);

    virtual VectorFiniteVolumeField compute(const std::vector<ImmersedBoundaryObject> &ibObjs);

protected:

    virtual void computeInterfaceNormals(const std::vector<ImmersedBoundaryObject> &ibObjs);

};

#endif

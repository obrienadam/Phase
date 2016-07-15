#ifndef IMMERSED_BOUNDARY_CELESTE_H
#define IMMERSED_BOUNDARY_CELESTE_H

#include "Celeste.h"

class ImmersedBoundaryObject;

class ImmersedBoundaryCeleste : public Celeste
{
public:

    ImmersedBoundaryCeleste(const Input& input,
                            Solver &solver);

    virtual VectorFiniteVolumeField compute(const std::vector<ImmersedBoundaryObject> &ibObjs);

protected:

    virtual void computeCurvature(const std::vector<ImmersedBoundaryObject> &ibObjs);
};

#endif

#ifndef IMMERSED_BOUNDARY_CELESTE_H
#define IMMERSED_BOUNDARY_CELESTE_H

#include "Celeste.h"

class ImmersedBoundaryObject;

class ImmersedBoundaryCeleste : public Celeste
{
public:

    ImmersedBoundaryCeleste(const Input& input,
                            const ScalarFiniteVolumeField& gamma,
                            const VectorFiniteVolumeField &u,
                            std::map<std::string, ScalarFiniteVolumeField>& scalarFields,
                            std::map<std::string, VectorFiniteVolumeField>& vectorFields);

    virtual VectorFiniteVolumeField compute(const std::vector<ImmersedBoundaryObject> &ibObjs);

protected:

    virtual void computeCurvature(const std::vector<ImmersedBoundaryObject> &ibObjs);
};

#endif

#ifndef CELESTE_H
#define CELESTE_H

#include "ContinuumSurfaceForce.h"

class Celeste : public ContinuumSurfaceForce
{
public:

    Celeste(const Input& input,
            const ScalarFiniteVolumeField& gamma,
            const VectorFiniteVolumeField &u,
            std::map<std::string, ScalarFiniteVolumeField>& scalarFields,
            std::map<std::string, VectorFiniteVolumeField>& vectorFields);

    virtual VectorFiniteVolumeField compute();

protected:

    virtual void computeGradGammaTilde();
    virtual void computeInterfaceNormals();
    virtual void computeCurvature();
    void weightCurvatures();

    ScalarFiniteVolumeField wGamma_;

};

#endif

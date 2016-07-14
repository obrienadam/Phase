#ifndef CELESTE_H
#define CELESTE_H

#include "ContinuumSurfaceForce.h"
#include "Matrix.h"

class Celeste : public ContinuumSurfaceForce
{
public:

    Celeste(const Input& input,
            const ScalarFiniteVolumeField& gamma,
            const VectorFiniteVolumeField &u,
            std::map<std::string, ScalarFiniteVolumeField>& scalarFields,
            std::map<std::string, VectorFiniteVolumeField>& vectorFields);

    virtual VectorFiniteVolumeField compute();
    void constructMatrices();

protected:

    virtual void computeGradGammaTilde();
    virtual void computeInterfaceNormals();
    virtual void computeCurvature();
    void weightCurvatures();

    std::vector<Matrix> matrices_;

    ScalarFiniteVolumeField wGamma_;

};

#endif

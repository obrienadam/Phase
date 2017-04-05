#ifndef CELESTE_H
#define CELESTE_H

#include "ContinuumSurfaceForce.h"
#include "Matrix.h"

class Celeste : public ContinuumSurfaceForce
{
public:

    Celeste(const Input& input,
            Solver &solver);

    virtual VectorFiniteVolumeField compute();

    void constructMatrices();

protected:

    void constructGammaTildeMatrices();
    void constructKappaMatrices();

    virtual void computeGradGammaTilde();
    virtual void computeInterfaceNormals();
    virtual void computeCurvature();
    void weightCurvatures();

    std::vector<Matrix> kappaMatrices_, gradGammaTildeMatrices_;
    ScalarFiniteVolumeField &w_;

};

#endif

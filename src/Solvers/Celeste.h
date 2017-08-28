#ifndef CELESTE_H
#define CELESTE_H

#include "SurfaceTensionForce.h"
#include "Matrix.h"

class Celeste : public SurfaceTensionForce
{
public:

    Celeste(const Input &input,
            const ImmersedBoundary &ib,
            const ScalarFiniteVolumeField &gamma,
            const ScalarFiniteVolumeField &rho,
            const ScalarFiniteVolumeField &mu,
            const VectorFiniteVolumeField &u,
            const ScalarGradient &gradGamma);

    void computeFaces();

    void computeFaces(const ImmersedBoundary& ib);

    void compute();

    void compute(const ImmersedBoundary& ib);

    void constructMatrices();

protected:

    virtual void computeGradGammaTilde();

    virtual void computeCurvature();

    void computeCurvature(const ImmersedBoundary& ib);

    Matrix leastSquaresMatrix(const Cell& cell, bool weighted = false);

    Matrix leastSquaresMatrix(const ImmersedBoundaryObject& ibObj, const Cell& cell, bool weighted = false);

    void constructGammaTildeMatrices();

    void constructKappaMatrices();

    std::vector<bool> modifiedStencil_;
    std::vector<Matrix> kappaMatrices_, gradGammaTildeMatrices_;
};

#endif

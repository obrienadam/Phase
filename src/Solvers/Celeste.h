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

    void compute();

    void compute(const ImmersedBoundary& ib);

    void constructMatrices();

    void constructMatrices(const ImmersedBoundary& ib);

protected:

    virtual void computeGradGammaTilde();

    virtual void computeCurvature();

    void computeCurvature(const ImmersedBoundary& ib);

    Matrix leastSquaresMatrix(const Cell& cell, bool weighted = false);

    Matrix leastSquaresMatrix(const ImmersedBoundary& ib, const Cell& cell, bool weighted = false);

    void constructGammaTildeMatrices();

    void constructKappaMatrices();

    std::vector<Matrix> kappaMatrices_, gradGammaTildeMatrices_;
};

#endif

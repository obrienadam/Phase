#ifndef MULTIPHASE_H
#define MULTIPHASE_H

#include "Simple.h"

class Multiphase : public Simple
{
public:
    Multiphase(const FiniteVolumeGrid2D& grid, const Input& input);

    ScalarFiniteVolumeField &gamma;

private:

    void computeRho();
    void computeMu();

    Scalar rho1_, rho2_, mu1_, mu2_, sigma_;
};

#endif

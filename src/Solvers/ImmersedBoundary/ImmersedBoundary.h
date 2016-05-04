#ifndef IMMERSED_BOUNDARY_H
#define IMMERSED_BOUNDARY_H

#include "Multiphase.h"
#include "ImmersedBoundaryObject.h"

class ImmersedBoundary : public Multiphase
{
public:

    enum {SOLID = 1, FLUID = 2, IB = 3};

    ImmersedBoundary(const FiniteVolumeGrid2D& grid, const Input& input);

    Scalar solve(Scalar timeStep, Scalar prevTimeStep);

    ScalarFiniteVolumeField &cellStatus_;

protected:

    void setCellStatus();

    Scalar solveUEqn(Scalar timeStep, Scalar prevTimeStep);
    Scalar solvePCorrEqn();
    Scalar solveGammaEqn(Scalar timeStep, Scalar prevTimeStep);

    ImmersedBoundaryObject ibObj_;
};

#endif

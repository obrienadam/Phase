#ifndef SIMPLE_H
#define SIMPLE_H

#include "Solver.h"
#include "Input.h"
#include "FiniteVolumeGrid2D.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"
#include "Equation.h"

class Simple : public Solver
{
public:

    Simple(const FiniteVolumeGrid2D& grid, const Input& input);

    virtual Scalar solve(Scalar timeStep) {}

    VectorFiniteVolumeField u, gradP;
    ScalarFiniteVolumeField p, pCorr;

protected:

    Equation<VectorFiniteVolumeField> uEqn_;
    Equation<ScalarFiniteVolumeField> pCorrEqn_;
};

#endif

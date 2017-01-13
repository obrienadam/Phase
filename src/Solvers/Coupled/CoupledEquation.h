#ifndef COUPLED_EQUATION_H
#define COUPLED_EQUATION_H

/**************************************************************
 *
 *
 * This class is deprecated and not currently functioning. May
 * be supported again in the future.
 *
 *
 **************************************************************/

#include "EigenSparseMatrixSolver.h"
#include "Vector.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"
/*
class CoupledEquation
{
public:

    CoupledEquation(const Input& input, const ScalarFiniteVolumeField& rho, const ScalarFiniteVolumeField& mu, VectorFiniteVolumeField& u, ScalarFiniteVolumeField& p);

    Scalar solve(Scalar timeStep);

protected:

    void assembleMomentumEquation(Scalar timeStep);
    void assembleContinuityEquation();
    void rhieChowInterpolation();

    size_t nActiveCells_, nVars_;

    const ScalarFiniteVolumeField &rho_, &mu_;
    VectorFiniteVolumeField& u_, gradP_;
    ScalarFiniteVolumeField& p_, d_;

    //const Vector2D& g_;

    std::vector<Triplet> triplets_;
    SparseMatrix spMat_;
    BiCGSTABDiagonalPreconditioner solver_;
    SparseVector rhs_, phi_;
};
*/
#endif

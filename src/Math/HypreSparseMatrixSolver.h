#ifndef HYPRE_SPARSE_MATRIX_SOLVER_H
#define HYPRE_SPARSE_MATRIX_SOLVER_H

#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_krylov.h>

#include "SparseMatrixSolver.h"
#include "Communicator.h"

class HypreSparseMatrixSolver : public SparseMatrixSolver
{
public:

    enum PreconditionerType
    {
        EUCLID, BOOMER_AMG
    };

    HypreSparseMatrixSolver(const Communicator &comm, PreconditionerType preconType = EUCLID);

    ~HypreSparseMatrixSolver();

    void setRank(int rank);

    void set(const CoefficientList &eqn);

    void setGuess(const Vector &x0);

    void setRhs(const Vector &rhs);

    Scalar solve();

    void mapSolution(ScalarFiniteVolumeField &field);

    void mapSolution(VectorFiniteVolumeField &field);

    void setMaxIters(int maxIters);

    void setToler(Scalar toler);

    void setDropToler(Scalar toler);

    void setFillFactor(int fill);

    int nIters() const
    { return nIters_; }

    Scalar error() const
    { return toler_; }

    bool supportsMPI() const
    { return true; }

    void printStatus(const std::string &msg) const;

    void init();

private:

    void initializeSolver();

    void deinitializeSolver();

    void initializeIJObjects(Size rank);

    void deinitializeIJObjects();

    PreconditionerType preconType_;

    HYPRE_Int nIters_, maxIters_, fill_;
    HYPRE_Complex minToler_, toler_, dropToler_;

    HYPRE_Int iLower_, iUpper_;
    std::vector<Size> localSizes_;
    std::vector<HYPRE_Int> globalInds_;

    const Communicator &comm_;

    HYPRE_IJMatrix ijMatrix_;
    HYPRE_IJVector b_, x_;

    HYPRE_Solver precon_;
    HYPRE_Solver solver_;
};

#endif

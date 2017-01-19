#ifndef HYPRE_SPARSE_MATRIX_SOLVER_H
#define HYPRE_SPARSE_MATRIX_SOLVER_H

#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>

#include "SparseMatrixSolver.h"
#include "Communicator.h"

class HypreSparseMatrixSolver : public SparseMatrixSolver
{
public:
    HypreSparseMatrixSolver(const Communicator &comm);
    ~HypreSparseMatrixSolver();

    void setRank(int rank);
    void set(const CoefficientList& eqn);
    void setRhs(const Vector& rhs);

    Scalar solve();

    void mapSolution(ScalarFiniteVolumeField& field);
    void mapSolution(VectorFiniteVolumeField& field);

    void setMaxIters(int maxIters) { maxIters_ = maxIters; }
    void setToler(Scalar toler) { minToler_ = toler; }
    void setFillFactor(int fill) { fill_ = fill; }

    int nIters() const { return nIters_; }
    Scalar error() const { return toler_; }

private:

    int nIters_, maxIters_, fill_;
    Scalar minToler_, toler_;

    int iLower_, iUpper_, jLower_, jUpper_;
    std::vector<Size> localSizes_;
    std::vector<int> globalInds_;

    const Communicator& comm_;

    HYPRE_IJMatrix ijMatrix_;
    HYPRE_ParCSRMatrix csrMatrix_;

    HYPRE_IJVector b_, x_;
    HYPRE_ParVector bPar_, xPar_;

    HYPRE_Solver solver_;
};

#endif

#ifndef PETSC_SPARSE_MATRIX_SOLVER_H
#define PETSC_SPARSE_MATRIX_SOLVER_H

#include <petsc.h>

#include "SparseMatrixSolver.h"
#include "Communicator.h"

class PetscSparseMatrixSolver : public SparseMatrixSolver
{
public:

    PetscSparseMatrixSolver(const Communicator &comm, const std::string &preconditioner);

    PetscSparseMatrixSolver(const PetscSparseMatrixSolver &rhs) = delete;

    ~PetscSparseMatrixSolver();

    void setRank(int rank);

    void set(const CoefficientList &eqn);

    void setGuess(const Vector &x0);

    void setRhs(const Vector &rhs);

    Scalar solve();

    void mapSolution(ScalarFiniteVolumeField &field);

    void mapSolution(VectorFiniteVolumeField &field);

    void setPreconditioner(const std::string &preconditioner);

    void setMaxIters(int maxIters);

    void setToler(Scalar toler);

    void setDropToler(Scalar toler);

    void setFillFactor(int fill);

    int nIters() const;

    Scalar error() const;

    bool supportsMPI() const;

    void printStatus(const std::string &msg) const;

private:

    static int solversActive_;

    const Communicator &comm_;
    PetscInt maxIters_ = 500, iLower_ = 0, iUpper_ = 0, error_;
    PetscReal rTol_ = 1e-10, absTol_ = 1e-10;

    Mat A_;
    Vec x_, b_;

    PC precon_;
    KSP solver_;
};

#endif

#ifndef SPARSE_MATRIX_SOLVER
#define SPARSE_MATRIX_SOLVER

#include "Vector.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

class SparseMatrixSolver
{
public:

    typedef std::pair<Index, Scalar> Entry;
    typedef std::vector<Entry> Row;
    typedef std::vector<Row> CoefficientList;

    virtual void setRank(int rank) = 0;
    virtual void set(const CoefficientList& eqn) = 0;
    virtual void setRhs(const Vector& rhs) = 0;

    virtual Scalar solve() = 0;
    virtual Scalar solve(const Vector& x0);

    virtual void mapSolution(ScalarFiniteVolumeField& field) = 0;
    virtual void mapSolution(VectorFiniteVolumeField& field) = 0;

    virtual void setPreconditioner(const std::string& preconditioner) {}
    virtual void setMaxIters(int maxIters) = 0;
    virtual void setToler(Scalar toler) = 0;
    virtual void setDropToler(Scalar toler) = 0;
    virtual void setFillFactor(int fill) = 0;
    virtual void setMaxPreconditionerUses(int maxPreconditionerUses);

    virtual int nIters() const = 0;
    virtual Scalar error() const = 0;

    virtual bool supportsMPI() const = 0;
    virtual void printStatus(const std::string& msg) const;

protected:
    int nPreconUses_ = 1, maxPreconUses_ = 1;
};

#ifdef PHASE_USE_EIGEN
#include "EigenSparseMatrixSolver.h"
#endif

#ifdef PHASE_USE_HYPRE
#include "HypreSparseMatrixSolver.h"
#endif

#ifdef PHASE_USE_PETSC
#include "PetscSparseMatrixSolver.h"
#endif

#endif

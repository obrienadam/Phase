#ifndef TRILINOS_BELOS_SPARSE_MATRIX_SOLVER_H
#define TRILINOS_BELOS_SPARSE_MATRIX_SOLVER_H

#include <Teuchos_DefaultMpiComm.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <BelosBiCGStabSolMgr.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Ifpack2_Factory.hpp>

#include "SparseMatrixSolver.h"

class TrilinosBelosSparseMatrixSolver : public SparseMatrixSolver
{
public:
    TrilinosBelosSparseMatrixSolver(const Communicator &comm);

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

    int nIters() const;

    Scalar error() const;

    bool supportsMPI() const
    { return true; }

    void printStatus(const std::string &msg) const;

private:

    typedef Teuchos::MpiComm<Index> TeuchosComm;
    typedef Tpetra::Map<Index, Index> TpetraMap;
    typedef Tpetra::RowMatrix<Scalar, Index, Index> TpetraRowMatrix;
    typedef Tpetra::CrsMatrix<Scalar, Index, Index> TpetraCrsMatrix;
    typedef Tpetra::Vector<Scalar, Index, Index> TpetraVector;
    typedef Tpetra::MultiVector<Scalar, Index, Index> TpetraMultiVector;
    typedef Tpetra::Operator<Scalar, Index, Index> Operator;
    typedef Belos::LinearProblem<Scalar, TpetraMultiVector, Operator> LinearProblem;
    typedef Belos::BiCGStabSolMgr<Scalar, TpetraMultiVector, Operator> Solver;
    typedef Ifpack2::RILUK<TpetraRowMatrix> Preconditioner;

    //- Communication objects
    const Communicator &comm_;
    Teuchos::RCP<TeuchosComm> Tcomm_;
    Teuchos::RCP<TpetraMap> map_;

    //- Parameters
    Teuchos::RCP<Teuchos::ParameterList> belosParameters_, ifpackParameters_;

    //- Matrix data structures
    Teuchos::RCP<TpetraCrsMatrix> mat_;
    Teuchos::RCP<TpetraVector> x_, b_;

    //- Solver data structures
    Teuchos::RCP<LinearProblem> linearProblem_;
    Teuchos::RCP<Solver> solver_;
    Teuchos::RCP<Preconditioner> precon_;
};

#endif

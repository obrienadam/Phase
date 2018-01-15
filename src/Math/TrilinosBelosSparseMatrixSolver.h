#ifndef TRILINOS_BELOS_SPARSE_MATRIX_SOLVER_H
#define TRILINOS_BELOS_SPARSE_MATRIX_SOLVER_H

#include <Teuchos_DefaultMpiComm.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosSolverManager.hpp>
#include <Ifpack2_Preconditioner.hpp>

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

    void setup(const boost::property_tree::ptree& parameters);

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
    typedef Tpetra::Operator<Scalar, Index, Index> TpetraOperator;
    typedef Belos::LinearProblem<Scalar, TpetraMultiVector, TpetraOperator> LinearProblem;
    typedef Belos::SolverManager<Scalar, TpetraMultiVector, TpetraOperator> Solver;
    typedef Ifpack2::Preconditioner<Scalar, Index, Index> Preconditioner;

    //- Communication objects
    const Communicator &comm_;
    Teuchos::RCP<TeuchosComm> Tcomm_;
    Teuchos::RCP<const TpetraMap> map_;

    //- Types
    std::string precType_;

    //- Parameters
    Teuchos::RCP<Teuchos::ParameterList> belosParams_, ifpackParams_;

    //- Factories

    //- Matrix data structures
    Teuchos::RCP<TpetraCrsMatrix> mat_;
    Teuchos::RCP<TpetraVector> x_, b_;

    //- Solver data structures
    Teuchos::RCP<LinearProblem> linearProblem_;
    Teuchos::RCP<Solver> solver_;
    Teuchos::RCP<Preconditioner> precon_;
};

#endif

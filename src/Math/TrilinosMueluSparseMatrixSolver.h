#ifndef TRILINOS_MUELU_SPARSE_MATRIX_SOLVER_H
#define TRILINOS_MUELU_SPARSE_MATRIX_SOLVER_H

#include <Teuchos_DefaultMpiComm.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <BelosTpetraAdapter.hpp>

#include "SparseMatrixSolver.h"

class TrilinosMueluSparseMatrixSolver: public SparseMatrixSolver
{
public:

    TrilinosMueluSparseMatrixSolver(const Communicator& comm);

    void setRank(int rank);

    void set(const CoefficientList &eqn);

    void setGuess(const Vector &x0);

    void setRhs(const Vector &rhs);

    Scalar solve();

    Scalar solve(const Vector &x0);

    void mapSolution(ScalarFiniteVolumeField &field);

    void mapSolution(VectorFiniteVolumeField &field);

    void setup(const boost::property_tree::ptree& parameters);

    int nIters() const;

    Scalar error() const;

    bool supportsMPI() const
    { return true; }

private:

    typedef Teuchos::MpiComm<Index> TeuchosComm;
    typedef Tpetra::Map<Index, Index> TpetraMap;
    typedef Tpetra::CrsMatrix<Scalar, Index, Index> TpetraCrsMatrix;
    typedef Tpetra::Vector<Scalar, Index, Index> TpetraVector;
    typedef Tpetra::MultiVector<Scalar, Index, Index> TpetraMultiVector;
    typedef Tpetra::Operator<Scalar, Index, Index> TpetraOperator;
    typedef Belos::LinearProblem<Scalar, TpetraMultiVector, TpetraOperator> LinearProblem;
    typedef Belos::SolverManager<Scalar, TpetraMultiVector, TpetraOperator> Solver;
    typedef MueLu::TpetraOperator<Scalar, Index, Index> MueLuTpetraOperator;

    const Communicator& comm_;
    Teuchos::RCP<TeuchosComm> Tcomm_;
    Teuchos::RCP<TpetraMap> map_;

    Teuchos::RCP<Teuchos::ParameterList> mueluParams_, belosParams_;
    Teuchos::RCP<MueLuTpetraOperator> precon_;

    Teuchos::RCP<TpetraVector> x_, b_;
    Teuchos::RCP<TpetraCrsMatrix> mat_;

    Teuchos::RCP<LinearProblem> linearProblem_;
    Teuchos::RCP<Solver> solver_;
};

#endif
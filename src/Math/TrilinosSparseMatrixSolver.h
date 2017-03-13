#ifndef TRILINOS_SPARSE_MATRIX_SOLVER_H
#define TRILINOS_SPARSE_MATRIX_SOLVER_H

#include <Teuchos_DefaultMpiComm.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Ifpack2_RILUK.hpp>
#include <Ifpack2_Diagonal.hpp>

#include "SparseMatrixSolver.h"

class TrilinosSparseMatrixSolver : public SparseMatrixSolver
{
public:
    TrilinosSparseMatrixSolver(const Communicator &comm,
                               const std::string &solver = "BiCGSTAB",
                               const std::string &preconType = "RILUK");

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

    void setMaxPreconditionerUses(int maxPreconditionerUses);

    int nIters() const;

    Scalar error() const;

    bool supportsMPI() const
    { return true; }

    void printStatus(const std::string &msg) const;

private:

    typedef const Teuchos::MpiComm<int> TeuchosComm;
    typedef const Tpetra::Map<Index, Index> TpetraMap;

    typedef Tpetra::CrsMatrix<Scalar, Index, Index> TpetraMatrix;
    typedef Tpetra::Vector<Scalar, Index, Index> TpetraVector;
    typedef Tpetra::MultiVector<Scalar, Index, Index> TpetraMultiVector;

    typedef Belos::LinearProblem<Scalar, TpetraMultiVector, Tpetra::Operator<Scalar, Index, Index>> LinearProblem;
    typedef Belos::SolverFactory<Scalar, TpetraMultiVector, Tpetra::Operator<Scalar, Index, Index>> SolverFactory;
    typedef Belos::SolverManager<Scalar, TpetraMultiVector, Tpetra::Operator<Scalar, Index, Index>> Solver;

    typedef Ifpack2::Preconditioner<Scalar, Index, Index> Preconditioner;

    const Communicator &comm_;
    Teuchos::RCP<TeuchosComm> Tcomm_;
    Teuchos::RCP<TpetraMap> map_;

    Teuchos::RCP<Teuchos::ParameterList> belosParameters_, ifpackParameters_;

    Teuchos::RCP<TpetraMatrix> mat_;
    Teuchos::RCP<TpetraVector> x_, b_;

    SolverFactory solverFactory_;
    Teuchos::RCP<LinearProblem> linearProblem_;
    std::string solverType_;
    Teuchos::RCP<Solver> solver_;

    std::string preconType_;
    //Ifpack2::Factory preconditionerFactory_;
    Teuchos::RCP<Preconditioner> preconditioner_;
};

#endif

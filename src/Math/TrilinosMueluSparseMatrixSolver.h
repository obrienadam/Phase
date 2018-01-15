#ifndef TRILINOS_MUELU_SPARSE_MATRIX_SOLVER_H
#define TRILINOS_MUELU_SPARSE_MATRIX_SOLVER_H

#include <Teuchos_DefaultMpiComm.hpp>
#include <Xpetra_TpetraRowMatrix.hpp>
#include <MueLu.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosSolverManager.hpp>

#include "SparseMatrixSolver.h"

class TrilinosMueluSparseMatrixSolver: public SparseMatrixSolver
{
public:

    TrilinosMueluSparseMatrixSolver(const Communicator& comm);

    TrilinosMueluSparseMatrixSolver(const Communicator& comm, const std::weak_ptr<const FiniteVolumeGrid2D>& grid);

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

    void printStatus(const std::string& msg) const;

private:

    typedef Teuchos::MpiComm<Index> TeuchosComm;

    typedef Tpetra::Map<Index, Index> TpetraMap;
    typedef Tpetra::Operator<Scalar, Index, Index> TpetraOperator;
    typedef Tpetra::CrsMatrix<Scalar, Index, Index> TpetraCrsMatrix;
    typedef Tpetra::Vector<Scalar, Index, Index> TpetraVector;
    typedef Tpetra::MultiVector<Scalar, Index, Index> TpetraMultiVector;
//
//    typedef Xpetra::Map<Index, Index> XpetraMap;
//    typedef Xpetra::CrsMatrix<Scalar, Index, Index> XpetraCrsMatrix;
//    typedef Xpetra::Vector<Scalar, Index, Index> XpetraVector;
//    typedef Xpetra::MultiVector<Scalar, Index, Index> XpetraMultiVector;


    typedef MueLu::Hierarchy<Scalar, Index, Index> MueLuHierarchy;
    typedef MueLu::HierarchyManager<Scalar, Index, Index> MueLuHierarchyManager;

    typedef Belos::LinearProblem<Scalar, TpetraMultiVector, TpetraOperator> LinearProblem;
    typedef Belos::SolverManager<Scalar, TpetraMultiVector, TpetraOperator> Solver;
    typedef MueLu::TpetraOperator<Scalar, Index, Index> Preconditioner;

    const Communicator& comm_;

    Teuchos::RCP<TeuchosComm> Tcomm_;
    Teuchos::RCP<TpetraMap> map_;

    Teuchos::RCP<Teuchos::ParameterList> belosParams_, mueluParams_;
    Teuchos::RCP<MueLuHierarchyManager> mueLuFactory_;
    Teuchos::RCP<MueLuHierarchy> h_;

    Teuchos::RCP<TpetraMultiVector> x_, b_, coords_;
    Teuchos::RCP<TpetraCrsMatrix> mat_;

    Teuchos::RCP<LinearProblem> linearProblem_;
    Teuchos::RCP<Solver> solver_;
    Teuchos::RCP<Preconditioner> precon_;

    std::weak_ptr<const FiniteVolumeGrid2D> grid_;
};

#endif
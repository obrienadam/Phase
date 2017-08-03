#ifndef TRILINOS_MUELU_SPARSE_MATRIX_SOLVER_H
#define TRILINOS_MUELU_SPARSE_MATRIX_SOLVER_H

#include <Teuchos_DefaultMpiComm.hpp>
#include <MueLu.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "SparseMatrixSolver.h"

class TrilinosMueluSparseMatrixSolver: public SparseMatrixSolver
{
public:

    TrilinosMueluSparseMatrixSolver(const Communicator& comm, const std::string& xmlFileName);

    void setRank(int rank);

    void set(const CoefficientList &eqn);

    void setGuess(const Vector &x0);

    void setRhs(const Vector &rhs);

    Scalar solve();

    Scalar solve(const Vector &x0);

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

private:

    typedef Teuchos::MpiComm<Index> TeuchosComm;
    typedef Tpetra::Map<Index, Index> TpetraMap;
    typedef Tpetra::RowMatrix<Scalar, Index, Index> TpetraRowMatrix;
    typedef Tpetra::CrsMatrix<Scalar, Index, Index> TpetraCrsMatrix;
    typedef Tpetra::Vector<Scalar, Index, Index> TpetraVector;
    typedef Tpetra::MultiVector<Scalar, Index, Index> TpetraMultiVector;
    typedef Tpetra::Operator<Scalar, Index, Index> Operator;
    typedef Belos::LinearProblem<Scalar, TpetraMultiVector, Operator> LinearProblem;

    const Communicator& comm_;
    Teuchos::RCP<TeuchosComm> Tcomm_;
    Teuchos::RCP<TpetraMap> map_;

    Teuchos::RCP<Teuchos::ParameterList> mueluParameters_;

};

#endif
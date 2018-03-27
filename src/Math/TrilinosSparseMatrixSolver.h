#ifndef TRILINOS_SPARSE_MATRIX_SOLVER_H
#define TRILINOS_SPARSE_MATRIX_SOLVER_H

#include <Tpetra_CrsMatrix.hpp>

#include "System/Communicator.h"

#include "SparseMatrixSolver.h"

class TrilinosSparseMatrixSolver : public SparseMatrixSolver
{
public:

    TrilinosSparseMatrixSolver(const Communicator &comm);

    virtual void setRank(int rank);

    virtual void set(const CoefficientList &eqn);

    virtual void setGuess(const Vector &x0);

    virtual void setRhs(const Vector &rhs);

    Scalar x(Index idx) const;

    virtual bool supportsMPI() const
    { return true; }

    virtual void printStatus(const std::string &msg) const;

protected:

    typedef Teuchos::MpiComm<Index> TeuchosComm;
    typedef Tpetra::Map<Index, Index> TpetraMap;
    typedef Tpetra::Operator<Scalar, Index, Index> TpetraOperator;
    typedef Tpetra::CrsMatrix<Scalar, Index, Index> TpetraCrsMatrix;
    typedef Tpetra::MultiVector<Scalar, Index, Index> TpetraMultiVector;

    const Communicator &comm_;

    Teuchos::RCP<TeuchosComm> Tcomm_;

    Teuchos::RCP<TpetraMap> map_;

    Teuchos::RCP<TpetraMultiVector> x_, b_;

    Teuchos::RCP<TpetraCrsMatrix> mat_;
};

#endif
#ifndef TRILINOS_SPARSE_MATRIX_SOLVER_H
#define TRILINOS_SPARSE_MATRIX_SOLVER_H

#include <Tpetra_CrsMatrix.hpp>

#include "System/Communicator.h"

#include "SparseMatrixSolver.h"

class TrilinosSparseMatrixSolver : public SparseMatrixSolver
{
public:

    typedef Teuchos::MpiComm<Index> TeuchosComm;
    typedef Tpetra::Map<Index, Index> TpetraMap;
    typedef Tpetra::Operator<Scalar, Index, Index> TpetraOperator;
    typedef Tpetra::CrsMatrix<Scalar, Index, Index> TpetraCrsMatrix;
    typedef Tpetra::MultiVector<Scalar, Index, Index> TpetraMultiVector;

    TrilinosSparseMatrixSolver(const Communicator &comm,
                               Tpetra::ProfileType pftype = Tpetra::StaticProfile);

    TrilinosSparseMatrixSolver(const Communicator &comm,  const Teuchos::RCP<TpetraCrsMatrix> &mat);

    virtual void setRank(int rank);

    virtual void setRank(int rowRank, int colRank);

    virtual void set(const CoefficientList &eqn) override;

    virtual void set(const std::vector<Index> &rowPtr, const std::vector<Index> &colInds, const std::vector<Scalar> &vals) override;

    virtual void set(const std::vector<SparseEntry> &entries) override;

    virtual void setGuess(const Vector &x0);

    virtual void setRhs(const Vector &rhs);

    virtual Scalar solveLeastSquares();

    Scalar x(Index idx) const
    { return xData_[idx]; }

    const Communicator &comm() const
    { return comm_; }

    virtual bool supportsMPI() const
    { return true; }

    virtual void printStatus(const std::string &msg) const;

    const Teuchos::RCP<TpetraCrsMatrix> &mat() const
    { return mat_; }

protected:

    const Communicator &comm_;

    Teuchos::RCP<const TeuchosComm> Tcomm_;

    Tpetra::ProfileType pftype_;

    Teuchos::RCP<TpetraMultiVector> x_, b_;

    Teuchos::RCP<const TpetraMap> rangeMap_, domainMap_;

    Teuchos::RCP<TpetraCrsMatrix> mat_;

    Teuchos::ArrayRCP<const Scalar> xData_;
};

std::shared_ptr<TrilinosSparseMatrixSolver> multiply(const TrilinosSparseMatrixSolver &A, const TrilinosSparseMatrixSolver &B, bool transA = false, bool transB = false);

#endif

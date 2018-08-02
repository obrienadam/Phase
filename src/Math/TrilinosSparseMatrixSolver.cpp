#include <TpetraExt_MatrixMatrix.hpp>

#include "TrilinosSparseMatrixSolver.h"

TrilinosSparseMatrixSolver::TrilinosSparseMatrixSolver(const Communicator &comm, Tpetra::ProfileType pftype)
    :
      comm_(comm),
      pftype_(pftype)
{
    Tcomm_ = rcp(new TeuchosComm(comm.communicator()));
}

TrilinosSparseMatrixSolver::TrilinosSparseMatrixSolver(const Communicator &comm,
                                                       const Teuchos::RCP<TrilinosSparseMatrixSolver::TpetraCrsMatrix> &mat)
    :
      comm_(comm),
      pftype_(mat->getProfileType()),
      mat_(mat)
{
    Tcomm_ = Teuchos::rcp_dynamic_cast<const TeuchosComm>(mat_->getComm(), true);
    x_ = rcp(new TpetraMultiVector(mat->getDomainMap(), 1, true));
    b_ = rcp(new TpetraMultiVector(mat->getRangeMap(), 1, true));
    xData_ = x_->getData(0);
}

void TrilinosSparseMatrixSolver::setRank(int rank)
{
    using namespace Teuchos;

    auto map = rcp(new TpetraMap(OrdinalTraits<Tpetra::global_size_t>::invalid(), rank, 0, Tcomm_));

    if (mat_.is_null() || !mat_->getRangeMap()->isSameAs(*map) || !mat_->getDomainMap()->isSameAs(*map)) //- Check if a new map is needed
    {
        rangeMap_ = domainMap_ = map;
        x_ = rcp(new TpetraMultiVector(map, 1, true));
        b_ = rcp(new TpetraMultiVector(map, 1, true));
        xData_ = x_->getData(0);
    }

    mat_ = rcp(new TpetraCrsMatrix(map, 20, pftype_));
}

void TrilinosSparseMatrixSolver::setRank(int rowRank, int colRank)
{
    using namespace Teuchos;

    auto rangeMap = rcp(new TpetraMap(OrdinalTraits<Tpetra::global_size_t>::invalid(), rowRank, 0, Tcomm_));
    auto domainMap = rcp(new TpetraMap(OrdinalTraits<Tpetra::global_size_t>::invalid(), colRank, 0, Tcomm_));

    if (mat_.is_null() || !mat_->getRangeMap()->isSameAs(*rangeMap) || !mat_->getDomainMap()->isSameAs(*domainMap))
    {
        rangeMap_ = rangeMap;
        domainMap_ = domainMap;
        x_ = rcp(new TpetraMultiVector(domainMap, 1, true));
        b_ = rcp(new TpetraMultiVector(rangeMap, 1, true));
        xData_ = x_->getData(0);
    }

    mat_ = rcp(new TpetraCrsMatrix(rangeMap, 9, pftype_));
}

void TrilinosSparseMatrixSolver::set(const CoefficientList &eqn)
{
    using namespace Teuchos;

    mat_->resumeFill();
    mat_->setAllToScalar(0.);

    std::vector<Index> cols; //- profiling shows that these should be outside
    std::vector<Scalar> vals;

    Index minGlobalIndex = mat_->getRowMap()->getMinGlobalIndex();

    for (Index localRow = 0, nLocalRows = eqn.size(); localRow < nLocalRows; ++localRow)
    {
        cols.clear();
        vals.clear();

        for (const auto &entry: eqn[localRow])
        {
            cols.push_back(entry.first);
            vals.push_back(entry.second);
        }

        mat_->insertGlobalValues(localRow + minGlobalIndex, cols.size(), vals.data(), cols.data());
    }

    mat_->fillComplete(domainMap_, rangeMap_);
}

void TrilinosSparseMatrixSolver::setGuess(const Vector &x0)
{
    x_->getDataNonConst(0).assign(std::begin(x0.data()), std::end(x0.data()));
}

void TrilinosSparseMatrixSolver::setRhs(const Vector &rhs)
{
    b_->getDataNonConst(0).assign(std::begin(rhs.data()), std::end(rhs.data()));
}

Scalar TrilinosSparseMatrixSolver::solveLeastSquares()
{
    auto A = mat_;

    auto C = Teuchos::rcp(new TpetraCrsMatrix(A->getDomainMap(), 12));
    auto b = Teuchos::rcp(new TpetraMultiVector(A->getDomainMap(), 1, false));

    Tpetra::MatrixMatrix::Multiply(*A, true, *A, false, *C, true);

    A->apply(*b_, *b, Teuchos::TRANS);

    mat_ = C;
    b_ = b;
    //- x_ should already have the correct domain map

    domainMap_ = mat_->getDomainMap();
    rangeMap_ = mat_->getRangeMap();

    solve();
}

void TrilinosSparseMatrixSolver::printStatus(const std::string &msg) const
{
    comm_ << msg << " iterations = " << nIters() << ", error = " << error() << ".\n";
}

//- External

std::shared_ptr<TrilinosSparseMatrixSolver> multiply(const TrilinosSparseMatrixSolver &A, const TrilinosSparseMatrixSolver &B, bool transA, bool transB)
{

}

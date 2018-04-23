#include "TrilinosSparseMatrixSolver.h"

TrilinosSparseMatrixSolver::TrilinosSparseMatrixSolver(const Communicator &comm)
        :
        comm_(comm)
{
    Tcomm_ = rcp(new TeuchosComm(comm.communicator()));
}

void TrilinosSparseMatrixSolver::setRank(int rank)
{
    using namespace Teuchos;

    auto map = rcp(new TpetraMap(OrdinalTraits<Tpetra::global_size_t>::invalid(), rank, 0, Tcomm_));

    if (map_.is_null() || !map_->isSameAs(*map)) //- Check if a new map is needed
    {
        map_ = map;
        x_ = rcp(new TpetraMultiVector(map_, 1, true));
        b_ = rcp(new TpetraMultiVector(map_, 1, true));
        xData_ = x_->getData(0);
    }

    mat_ = rcp(new TpetraCrsMatrix(map_, 9, Tpetra::StaticProfile));
}

void TrilinosSparseMatrixSolver::set(const CoefficientList &eqn)
{
    using namespace Teuchos;

    Index minGlobalIndex = map_->getMinGlobalIndex();

    mat_->resumeFill();
    mat_->setAllToScalar(0.);

    std::vector<Index> cols; //- profiling shows that these should be outside
    std::vector<Scalar> vals;

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

    mat_->fillComplete();
}

void TrilinosSparseMatrixSolver::setGuess(const Vector &x0)
{
    x_->getDataNonConst(0).assign(std::begin(x0.data()), std::end(x0.data()));
}

void TrilinosSparseMatrixSolver::setRhs(const Vector &rhs)
{
    b_->getDataNonConst(0).assign(std::begin(rhs.data()), std::end(rhs.data()));
}

void TrilinosSparseMatrixSolver::printStatus(const std::string &msg) const
{
    comm_ << msg << " iterations = " << nIters() << ", error = " << error() << ".\n";
}

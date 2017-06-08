#include <MueLu_CreateTpetraPreconditioner.hpp>

#include "TrilinosSparseMatrixSolver.h"

TrilinosSparseMatrixSolver::TrilinosSparseMatrixSolver(const Communicator &comm)
        :
        comm_(comm)
{
    Tcomm_ = rcp(new TeuchosComm(comm.communicator()));
    belosParameters_ = rcp(new Teuchos::ParameterList());
    ifpackParameters_ = rcp(new Teuchos::ParameterList());
    solver_ = rcp(new Solver(rcp(new LinearProblem()), belosParameters_));
}

void TrilinosSparseMatrixSolver::setRank(int rank)
{
    using namespace Teuchos;

    auto map = rcp(new TpetraMap(OrdinalTraits<Tpetra::global_size_t>::invalid(), rank, 0, Tcomm_));

    if(map_.is_null() || !map_->isSameAs(*map))
    {
        map_ = map;
        mat_ = rcp(new TpetraCrsMatrix(map_, 5, Tpetra::DynamicProfile));
        x_ = rcp(new TpetraVector(map_, true));
        b_ = rcp(new TpetraVector(map_, true));
        precon_ = rcp(new Preconditioner(rcp_static_cast<const TpetraRowMatrix>(mat_)));
        precon_->setParameters(*ifpackParameters_);
        nPreconUses_ = maxPreconUses_; //- Ensure preconditioner gets recomputed
        linearProblem_ = rcp(new LinearProblem(mat_, x_, b_));
        linearProblem_->setRightPrec(precon_);
        solver_->setProblem(linearProblem_);
    }
}

void TrilinosSparseMatrixSolver::set(const SparseMatrixSolver::CoefficientList &eqn)
{
    Index minGlobalIndex = map_->getMinGlobalIndex();

    mat_->resumeFill();
    mat_->setAllToScalar(0);
    for (Index localRow = 0, nLocalRows = eqn.size(); localRow < nLocalRows; ++localRow)
    {
        std::vector<Index> cols;
        std::vector<Scalar> vals;

        for (const auto &entry: eqn[localRow])
        {
            cols.push_back(entry.first);
            vals.push_back(entry.second);
        }

        if(mat_->getProfileType() == Tpetra::StaticProfile)
            mat_->replaceGlobalValues(localRow + minGlobalIndex, cols.size(), vals.data(), cols.data());
        else
            mat_->insertGlobalValues(localRow + minGlobalIndex, cols.size(), vals.data(), cols.data());
    }

    bool initPrecon = mat_->getProfileType() == Tpetra::DynamicProfile;

    mat_->fillComplete();

    if(initPrecon)
        precon_->initialize();
}

void TrilinosSparseMatrixSolver::setGuess(const Vector &x0)
{
    std::transform(x0.begin(), x0.end(), x_->getDataNonConst().begin(), [](Scalar val) { return val; });
}

void TrilinosSparseMatrixSolver::setRhs(const Vector &rhs)
{
    std::transform(rhs.begin(), rhs.end(), b_->getDataNonConst().begin(), [](Scalar val) { return val; });
}

Scalar TrilinosSparseMatrixSolver::solve()
{
    if(nPreconUses_++ == maxPreconUses_)
    {
        precon_->compute();
        nPreconUses_ = 1;
    }

    linearProblem_->setProblem();
    solver_->solve();

    return error();
}

void TrilinosSparseMatrixSolver::mapSolution(ScalarFiniteVolumeField &field)
{
    Teuchos::ArrayRCP<const Scalar> soln = x_->getData();
    for (const Cell &cell: field.grid.localActiveCells())
        field(cell) = soln[cell.index(0)];
}

void TrilinosSparseMatrixSolver::mapSolution(VectorFiniteVolumeField &field)
{
    Teuchos::ArrayRCP<const Scalar> soln = x_->getData();
    Index nActiveCells = field.grid.localActiveCells().size();

    for (const Cell &cell: field.grid.localActiveCells())
    {
        field(cell).x = soln[cell.index(0)];
        field(cell).y = soln[cell.index(0) + nActiveCells];
    }
}

void TrilinosSparseMatrixSolver::setMaxIters(int maxIters)
{
    belosParameters_->set("Maximum Iterations", maxIters);
    solver_->setParameters(belosParameters_);
}

void TrilinosSparseMatrixSolver::setToler(Scalar toler)
{
    belosParameters_->set("Convergence Tolerance", toler);
    solver_->setParameters(belosParameters_);
}

void TrilinosSparseMatrixSolver::setDropToler(Scalar toler)
{
    ifpackParameters_->set("fact: relax value", toler);
}

void TrilinosSparseMatrixSolver::setFillFactor(int fill)
{
    ifpackParameters_->set("fact: iluk level-of-fill", fill);
}

int TrilinosSparseMatrixSolver::nIters() const
{
    return solver_->getNumIters();
}

Scalar TrilinosSparseMatrixSolver::error() const
{
    return solver_->achievedTol();
}

void TrilinosSparseMatrixSolver::printStatus(const std::string &msg) const
{
    comm_.printf("%s %s iterations = %d, error = %lf.\n", msg.c_str(), "BiCGSTAB", nIters(), error());
}

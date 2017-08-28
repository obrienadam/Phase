#include <MueLu_CreateTpetraPreconditioner.hpp>
#include "TrilinosMueluSparseMatrixSolver.h"

TrilinosMueluSparseMatrixSolver::TrilinosMueluSparseMatrixSolver(const Communicator &comm)
        :
        comm_(comm)
{
    Tcomm_ = rcp(new TeuchosComm(comm.communicator()));
    belosParams_ = rcp(new Teuchos::ParameterList());
    mueluParams_ = rcp(new Teuchos::ParameterList());
}

void TrilinosMueluSparseMatrixSolver::setRank(int rank)
{
    using namespace Teuchos;

    auto map = rcp(new TpetraMap(OrdinalTraits<Tpetra::global_size_t>::invalid(), rank, 0, Tcomm_));

    if (map_.is_null() || !map_->isSameAs(*map)) //- Check if a new map is needed
    {
        map_ = map;
        mat_ = rcp(new TpetraCrsMatrix(map_, 5, Tpetra::StaticProfile));
        x_ = rcp(new TpetraVector(map_, true));
        b_ = rcp(new TpetraVector(map_, true));

        std::cout << mueluParams_.is_null() << std::endl;

        precon_ = MueLu::CreateTpetraPreconditioner(rcp_static_cast<TpetraOperator>(mat_), "mg.xml");
        nPreconUses_ = maxPreconUses_; //- Ensure preconditioner gets recomputed

        linearProblem_ = rcp(new LinearProblem(mat_, x_, b_));
        linearProblem_->setLeftPrec(precon_);

        solver_->setProblem(linearProblem_);
        solver_->setParameters(belosParams_);
    }
}

void TrilinosMueluSparseMatrixSolver::set(const CoefficientList &eqn)
{
    using namespace Teuchos;

    Index minGlobalIndex = map_->getMinGlobalIndex();

    bool newMat = !mat_->isFillComplete();

    mat_->resumeFill();
    mat_->setAllToScalar(0.);
    for (Index localRow = 0, nLocalRows = eqn.size(); localRow < nLocalRows; ++localRow)
    {
        std::vector<Index> cols;
        std::vector<Scalar> vals;

        for (const auto &entry: eqn[localRow])
        {
            cols.push_back(entry.first);
            vals.push_back(entry.second);
        }

        if (newMat)
            mat_->insertGlobalValues(localRow + minGlobalIndex, cols.size(), vals.data(), cols.data());
        else
            mat_->replaceGlobalValues(localRow + minGlobalIndex, cols.size(), vals.data(), cols.data());
    }

    mat_->fillComplete();
}

void TrilinosMueluSparseMatrixSolver::setGuess(const Vector &x0)
{
    std::transform(x0.begin(), x0.end(), x_->getDataNonConst().begin(), [](Scalar val) { return val; });
}

void TrilinosMueluSparseMatrixSolver::setRhs(const Vector &rhs)
{
    std::transform(rhs.begin(), rhs.end(), b_->getDataNonConst().begin(), [](Scalar val) { return val; });
}

Scalar TrilinosMueluSparseMatrixSolver::solve()
{
    comm_.printf("Belos: Performing Krylov iterations...\n");
    linearProblem_->setProblem();

    try
    {
        solver_->solve();
    }
    catch (const Belos::StatusTestError &e)
    {
        comm_.printf("Error detected! Setting solution vector to 0 and attempting to resolve...\n");
        comm_.barrier();
        x_->putScalar(0.);
        solve();
    }

    return error();
}

Scalar TrilinosMueluSparseMatrixSolver::solve(const Vector &x0)
{

}

void TrilinosMueluSparseMatrixSolver::mapSolution(ScalarFiniteVolumeField &field)
{
    Teuchos::ArrayRCP<const Scalar> soln = x_->getData();
    for (const Cell &cell: field.grid().localActiveCells())
        field(cell) = soln[cell.index(0)];
}

void TrilinosMueluSparseMatrixSolver::mapSolution(VectorFiniteVolumeField &field)
{
    Teuchos::ArrayRCP<const Scalar> soln = x_->getData();
    Index nActiveCells = field.grid().localActiveCells().size();

    for (const Cell &cell: field.grid().localActiveCells())
    {
        field(cell).x = soln[cell.index(0)];
        field(cell).y = soln[cell.index(0) + nActiveCells];
    }
}

void TrilinosMueluSparseMatrixSolver::setup(const boost::property_tree::ptree &parameters)
{
    Factory factory;

    belosParams_->set("Maximum Iterations", parameters.get<int>("maxIters", 500));
    belosParams_->set("Convergence Tolerance", parameters.get<Scalar>("tolerance", 1e-8));
    solver_ = factory.create(parameters.get<std::string>("solver", "BICGSTAB"), belosParams_);

    std::cout << mueluParams_ << std::endl;
    mueluParams_ = Teuchos::getParametersFromXmlFile(parameters.get<std::string>("parameterFile"));
    std::cout << mueluParams_ << std::endl;
}

int TrilinosMueluSparseMatrixSolver::nIters() const
{
    return solver_->getNumIters();
}

Scalar TrilinosMueluSparseMatrixSolver::error() const
{
    return solver_->achievedTol();
}
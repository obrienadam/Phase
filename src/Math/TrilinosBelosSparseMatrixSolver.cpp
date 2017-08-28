#include <BelosSolverFactory.hpp>
#include <Ifpack2_Factory.hpp>

#include "TrilinosBelosSparseMatrixSolver.h"

TrilinosBelosSparseMatrixSolver::TrilinosBelosSparseMatrixSolver(const Communicator &comm)
        :
        comm_(comm)
{
    Tcomm_ = rcp(new TeuchosComm(comm.communicator()));
    belosParams_ = rcp(new Teuchos::ParameterList());
    ifpackParams_ = rcp(new Teuchos::ParameterList());
    schwarzParams_ = rcp(new Teuchos::ParameterList());
}

void TrilinosBelosSparseMatrixSolver::setRank(int rank)
{
    using namespace Teuchos;

    auto map = rcp(new TpetraMap(OrdinalTraits<Tpetra::global_size_t>::invalid(), rank, 0, Tcomm_));

    if (map_.is_null() || !map_->isSameAs(*map)) //- Check if a new map is needed
    {
        map_ = map;
        mat_ = rcp(new TpetraCrsMatrix(map_, 5, Tpetra::DynamicProfile));
        x_ = rcp(new TpetraVector(map_, true));
        b_ = rcp(new TpetraVector(map_, true));

        Ifpack2::Factory factory;

        precon_ = factory.create("SCHWARZ", rcp_static_cast<const TpetraCrsMatrix>(mat_));
        schwarzParams_->set("schwarz: inner preconditioner parameters", *ifpackParams_);
        precon_->setParameters(*schwarzParams_);
        nPreconUses_ = maxPreconUses_; //- Ensure preconditioner gets recomputed

        linearProblem_ = rcp(new LinearProblem(mat_, x_, b_));
        linearProblem_->setRightPrec(precon_);
        solver_->setProblem(linearProblem_);
        solver_->setParameters(belosParams_);
    }
}

void TrilinosBelosSparseMatrixSolver::set(const SparseMatrixSolver::CoefficientList &eqn)
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

        if(newMat)
            mat_->insertGlobalValues(localRow + minGlobalIndex, cols.size(), vals.data(), cols.data());
        else
            mat_->replaceGlobalValues(localRow + minGlobalIndex, cols.size(), vals.data(), cols.data());
    }

    mat_->fillComplete();
}

void TrilinosBelosSparseMatrixSolver::setGuess(const Vector &x0)
{
    std::transform(x0.begin(), x0.end(), x_->getDataNonConst().begin(), [](Scalar val) { return val; });
}

void TrilinosBelosSparseMatrixSolver::setRhs(const Vector &rhs)
{
    std::transform(rhs.begin(), rhs.end(), b_->getDataNonConst().begin(), [](Scalar val) { return val; });
}

Scalar TrilinosBelosSparseMatrixSolver::solve()
{
    if (nPreconUses_++ >= maxPreconUses_)
    {
        comm_.printf("Ifpack2: Computing preconditioner...\n");
        precon_->initialize();
        precon_->compute();
        nPreconUses_ = 1;
    }

    comm_.printf("Belos: Performing BiCGSTAB iterations...\n");
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

void TrilinosBelosSparseMatrixSolver::mapSolution(ScalarFiniteVolumeField &field)
{
    Teuchos::ArrayRCP<const Scalar> soln = x_->getData();
    for (const Cell &cell: field.grid().localActiveCells())
        field(cell) = soln[cell.index(0)];
}

void TrilinosBelosSparseMatrixSolver::mapSolution(VectorFiniteVolumeField &field)
{
    Teuchos::ArrayRCP<const Scalar> soln = x_->getData();
    Index nActiveCells = field.grid().localActiveCells().size();

    for (const Cell &cell: field.grid().localActiveCells())
    {
        field(cell).x = soln[cell.index(0)];
        field(cell).y = soln[cell.index(0) + nActiveCells];
    }
}

void TrilinosBelosSparseMatrixSolver::setup(const boost::property_tree::ptree& parameters)
{
    typedef Belos::SolverFactory<Scalar, TpetraMultiVector, Operator> SolverFactory;

    SolverFactory factory;

    belosParams_->set("Maximum Iterations", parameters.get<int>("maxIters", 500));
    belosParams_->set("Convergence Tolerance", parameters.get<Scalar>("tolerance", 1e-8));
    solver_ = factory.create(parameters.get<std::string>("solver", "BICGSTAB"), belosParams_);

    ifpackParams_->set("fact: iluk level-of-fill", parameters.get<Scalar>("iluFill", 0.));
    ifpackParams_->set("fact: ilut level-of-fill", parameters.get<Scalar>("iluFill", 1.));
    schwarzParams_->set("schwarz: inner preconditioner name", parameters.get<std::string>("innerPreconditioner", "RILUK"));
    schwarzParams_->set("schwarz: num iterations", parameters.get<int>("schwarzIters", 1));
    schwarzParams_->set("schwarz: combine mode", parameters.get<std::string>("schwarzCombineMode", "ADD"));
    schwarzParams_->set("schwarz: overlap level", parameters.get<int>("schwarzOverlap", 0));
}

int TrilinosBelosSparseMatrixSolver::nIters() const
{
    return solver_->getNumIters();
}

Scalar TrilinosBelosSparseMatrixSolver::error() const
{
    return solver_->achievedTol();
}

void TrilinosBelosSparseMatrixSolver::printStatus(const std::string &msg) const
{
    comm_.printf("%s %s iterations = %d, error = %lf.\n", msg.c_str(), "Krylov", nIters(), error());
}

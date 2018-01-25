#include <Teuchos_XMLParameterListCoreHelpers.hpp>
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
    linearProblem_ = rcp(new LinearProblem());
}

void TrilinosBelosSparseMatrixSolver::setRank(int rank)
{
    using namespace Teuchos;

    auto map = rcp(new const TpetraMap(OrdinalTraits<Tpetra::global_size_t>::invalid(), rank, 0, Tcomm_));

    if (map_.is_null() || !map_->isSameAs(*map)) //- Check if a new map is needed
    {
        map_ = map;
        x_ = rcp(new TpetraVector(map_, true));
        b_ = rcp(new TpetraVector(map_, false));
    }

    mat_ = rcp(new TpetraCrsMatrix(map_, 9, Tpetra::StaticProfile));
    precon_ = Ifpack2::Factory().create(precType_, rcp_static_cast<const TpetraRowMatrix>(mat_));
    precon_->setParameters(*ifpackParams_);
    linearProblem_->setOperator(mat_);
    linearProblem_->setRightPrec(precon_);
}

void TrilinosBelosSparseMatrixSolver::set(const SparseMatrixSolver::CoefficientList &eqn)
{
    using namespace Teuchos;

    Index minGlobalIndex = map_->getMinGlobalIndex();

    for (Index localRow = 0, nLocalRows = eqn.size(); localRow < nLocalRows; ++localRow)
    {
        std::vector<Index> cols;
        std::vector<Scalar> vals;

        for (const auto &entry: eqn[localRow])
        {
            cols.push_back(entry.first);
            vals.push_back(entry.second);
        }

        mat_->insertGlobalValues(localRow + minGlobalIndex, cols.size(), vals.data(), cols.data());
    }

    mat_->fillComplete();
}

void TrilinosBelosSparseMatrixSolver::setGuess(const Vector &x0)
{
    x_->getDataNonConst().assign(std::begin(x0.data()), std::end(x0.data()));
}

void TrilinosBelosSparseMatrixSolver::setRhs(const Vector &rhs)
{
    b_->getDataNonConst().assign(std::begin(rhs.data()), std::end(rhs.data()));
}

Scalar TrilinosBelosSparseMatrixSolver::solve()
{
    comm_.printf("Ifpack2: Computing preconditioner...\n");
    precon_->initialize();
    precon_->compute();

    comm_.printf("Belos: Performing Krylov iterations...\n");
    linearProblem_->setProblem(x_, b_);
    solver_->solve();

    return error();
}

void TrilinosBelosSparseMatrixSolver::mapSolution(ScalarFiniteVolumeField &field)
{
    Teuchos::ArrayRCP<const Scalar> soln = x_->getData();
    for (const Cell &cell: field.grid()->localActiveCells())
        field(cell) = soln[field.indexMap()->local(cell, 0)];
}

void TrilinosBelosSparseMatrixSolver::mapSolution(VectorFiniteVolumeField &field)
{
    Teuchos::ArrayRCP<const Scalar> soln = x_->getData();

    for (const Cell &cell: field.grid()->localActiveCells())
    {
        Vector2D &vec = field(cell);
        vec.x = soln[field.indexMap()->local(cell, 0)];
        vec.y = soln[field.indexMap()->local(cell, 1)];
    }
}

void TrilinosBelosSparseMatrixSolver::setup(const boost::property_tree::ptree& parameters)
{
    typedef Belos::SolverFactory<Scalar, TpetraMultiVector, TpetraOperator> SolverFactory;

    std::string filename = parameters.get<std::string>("belosParamFile", "");

    if(filename.empty())
    {
        belosParams_->set("Maximum Iterations", parameters.get<int>("maxIters", 500));
        belosParams_->set("Convergence Tolerance", parameters.get<Scalar>("tolerance", 1e-8));
    }
    else
        belosParams_ = Teuchos::getParametersFromXmlFile("case/" + filename);

    solver_ = SolverFactory().create(parameters.get<std::string>("solver", "BICGSTAB"), belosParams_);
    solver_->setProblem(linearProblem_);

    precType_ = parameters.get<std::string>("preconditioner", "schwarz");
    filename = parameters.get<std::string>("ifpackParamFile", "");

    if(filename.empty())
    {
        precType_ = "schwarz";

        auto tmp = rcp(new Teuchos::ParameterList());
        tmp->set("fact: iluk level-of-fill", parameters.get<Scalar>("iluFill", 0.));
        tmp->set("fact: ilut level-of-fill", parameters.get<Scalar>("iluFill", 1.));

        ifpackParams_->set("schwarz: inner preconditioner name", parameters.get<std::string>("innerPreconditioner", "RILUK"));
        ifpackParams_->set("schwarz: num iterations", parameters.get<int>("schwarzIters", 1));
        ifpackParams_->set("schwarz: combine mode", parameters.get<std::string>("schwarzCombineMode", "ADD"));
        ifpackParams_->set("schwarz: overlap level", parameters.get<int>("schwarzOverlap", 0));
        ifpackParams_->set("schwarz: inner preconditioner parameters", *tmp);
    }
    else
        ifpackParams_ = Teuchos::getParametersFromXmlFile("case/" + filename);
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

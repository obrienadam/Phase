#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <BelosSolverFactory.hpp>
#include <Ifpack2_Factory.hpp>

#include "TrilinosBelosSparseMatrixSolver.h"

TrilinosBelosSparseMatrixSolver::TrilinosBelosSparseMatrixSolver(const Communicator &comm)
    :
      TrilinosSparseMatrixSolver(comm)
{
    typedef Belos::SolverFactory<Scalar, TpetraMultiVector, TpetraOperator> SolverFactory;

    belosParams_ = rcp(new Teuchos::ParameterList());
    ifpackParams_ = rcp(new Teuchos::ParameterList());
    linearProblem_ = rcp(new LinearProblem());
}

void TrilinosBelosSparseMatrixSolver::setRank(int rank)
{
    using namespace Teuchos;
    typedef Tpetra::RowMatrix<Scalar, Index, Index> TpetraRowMatrix;

    TrilinosSparseMatrixSolver::setRank(rank);

    precon_ = Ifpack2::Factory().create(precType_, rcp_static_cast<const TpetraRowMatrix>(mat_));
    precon_->setParameters(*ifpackParams_);
    linearProblem_->setOperator(mat_);
    linearProblem_->setRightPrec(precon_);
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

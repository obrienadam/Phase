#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <BelosSolverFactory.hpp>

#include "TrilinosMueluSparseMatrixSolver.h"

TrilinosMueluSparseMatrixSolver::TrilinosMueluSparseMatrixSolver(const Communicator &comm)
        :
        TrilinosSparseMatrixSolver(comm)
{
//    h_ = Teuchos::rcp(new MueLuHierarchy());
//    h_->setDefaultVerbLevel(Teuchos::VERB_HIGH);
//    h_->GetLevel()->setDefaultVerbLevel(Teuchos::VERB_HIGH);
//    h_->IsPreconditioner(false);
}

TrilinosMueluSparseMatrixSolver::TrilinosMueluSparseMatrixSolver(const Communicator &comm,
                                                                 const std::weak_ptr<const FiniteVolumeGrid2D> &grid)
        :
        TrilinosSparseMatrixSolver(comm),
        grid_(grid)
{

}

void TrilinosMueluSparseMatrixSolver::setRank(int rank)
{
    using namespace Teuchos;

    TrilinosSparseMatrixSolver::setRank(rank);

    linearProblem_->setOperator(mat_);

    //h_->GetLevel(0)->Set("A", MueLu::TpetraCrs_To_XpetraMatrix(mat_));

    auto grid = grid_.lock();

    if (grid)
    {
        if (coords_.is_null() || !coords_->getMap()->isSameAs(*map_))
            coords_ = rcp(new TpetraMultiVector(map_, 2));

        std::transform(grid->localActiveCells().begin(),
                       grid->localActiveCells().end(),
                       coords_->getDataNonConst(0).begin(),
                       [](const Cell &cell) { return cell.centroid().x; });


        std::transform(grid->localActiveCells().begin(),
                       grid->localActiveCells().end(),
                       coords_->getDataNonConst(1).begin(),
                       [](const Cell &cell) { return cell.centroid().y; });
        //h_->GetLevel(0)->Set("Coordinates", MueLu::TpetraMultiVector_To_XpetraMultiVector(coords_));
    }
    else
        coords_ = Teuchos::null;
}

Scalar TrilinosMueluSparseMatrixSolver::solve()
{
//    comm_.printf("MueLu: Performing multigrid iterations...\n");
//
//    if (mueLuFactory_.is_null())
//        h_->Setup();
//    else
//        mueLuFactory_->SetupHierarchy(*h_);
//
//
//    h_->Iterate(*MueLu::TpetraMultiVector_To_XpetraMultiVector(b_),
//                *MueLu::TpetraMultiVector_To_XpetraMultiVector(x_));
    precon_ = MueLu::CreateTpetraPreconditioner(
            Teuchos::rcp_static_cast<Tpetra::Operator<Scalar, Index, Index>>(mat_),
            *mueluParams_,
            coords_);

    linearProblem_->setProblem(x_, b_);
    linearProblem_->setLeftPrec(precon_);

    solver_->solve();

    return error();
}

void TrilinosMueluSparseMatrixSolver::setup(const boost::property_tree::ptree &parameters)
{
    using namespace std;
    using namespace Teuchos;

    typedef Belos::SolverFactory<Scalar, TpetraMultiVector, TpetraOperator> SolverFactory;

    std::string belosParamFile = parameters.get<std::string>("belosParamFile");
    std::string mueluParamFile = parameters.get<std::string>("mueluParamFile");
    std::string solverName = parameters.get<std::string>("solver", "GMRES");

    belosParams_ = Teuchos::getParametersFromXmlFile("case/" + belosParamFile);
    mueluParams_ = Teuchos::getParametersFromXmlFile("case/" + mueluParamFile);
    solver_ = SolverFactory().create(solverName, belosParams_);

    linearProblem_ = rcp(new LinearProblem());
    solver_->setProblem(linearProblem_);
}

int TrilinosMueluSparseMatrixSolver::nIters() const
{
    return solver_->getNumIters();
}

Scalar TrilinosMueluSparseMatrixSolver::error() const
{
    return solver_->achievedTol();
}

void TrilinosMueluSparseMatrixSolver::printStatus(const std::string &msg) const
{
    comm_.printf("%s %s iterations = %d, error = %lf.\n", msg.c_str(), "Krylov", nIters(), error());
}
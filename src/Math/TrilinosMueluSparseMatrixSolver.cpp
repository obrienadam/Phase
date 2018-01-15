#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <BelosSolverFactory.hpp>

#include "TrilinosMueluSparseMatrixSolver.h"

TrilinosMueluSparseMatrixSolver::TrilinosMueluSparseMatrixSolver(const Communicator &comm)
        :
        comm_(comm)
{
    Tcomm_ = rcp(new TeuchosComm(comm.communicator()));
//    h_ = Teuchos::rcp(new MueLuHierarchy());
//    h_->setDefaultVerbLevel(Teuchos::VERB_HIGH);
//    h_->GetLevel()->setDefaultVerbLevel(Teuchos::VERB_HIGH);
//    h_->IsPreconditioner(false);
}

TrilinosMueluSparseMatrixSolver::TrilinosMueluSparseMatrixSolver(const Communicator &comm,
                                                                 const std::weak_ptr<const FiniteVolumeGrid2D> &grid)
        :
        TrilinosMueluSparseMatrixSolver(comm)
{
    grid_ = grid;
}

void TrilinosMueluSparseMatrixSolver::setRank(int rank)
{
    using namespace Teuchos;

    auto map = rcp(new TpetraMap(OrdinalTraits<Tpetra::global_size_t>::invalid(), rank, 0, Tcomm_));

    if (map_.is_null() || !map_->isSameAs(*map)) //- Check if a new map is needed
    {
        map_ = map;
        x_ = rcp(new TpetraMultiVector(map_, 1, true));
        b_ = rcp(new TpetraMultiVector(map_, 1, true));
        coords_ = rcp(new TpetraMultiVector(map_, 2));
    }

    mat_ = rcp(new TpetraCrsMatrix(map_, 8, Tpetra::StaticProfile));

    linearProblem_->setOperator(mat_);

    //h_->GetLevel(0)->Set("A", MueLu::TpetraCrs_To_XpetraMatrix(mat_));

    auto grid = grid_.lock();

    if (grid)
    {
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
}

void TrilinosMueluSparseMatrixSolver::set(const CoefficientList &eqn)
{
    using namespace Teuchos;

    Index minGlobalIndex = map_->getMinGlobalIndex();

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

        mat_->insertGlobalValues(localRow + minGlobalIndex, cols.size(), vals.data(), cols.data());
    }

    mat_->fillComplete();
}

void TrilinosMueluSparseMatrixSolver::setGuess(const Vector &x0)
{
    x_->getDataNonConst(0).assign(x0.begin(), x0.end());
}

void TrilinosMueluSparseMatrixSolver::setRhs(const Vector &rhs)
{
    b_->getDataNonConst(0).assign(rhs.begin(), rhs.end());
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

Scalar TrilinosMueluSparseMatrixSolver::solve(const Vector &x0)
{
    setGuess(x0);
    solve();
}

void TrilinosMueluSparseMatrixSolver::mapSolution(ScalarFiniteVolumeField &field)
{
    Teuchos::ArrayRCP<const Scalar> soln = x_->getData(0);
    for (const Cell &cell: field.grid().localActiveCells())
        field(cell) = soln[field.indexMap()->local(cell, 0)];
}

void TrilinosMueluSparseMatrixSolver::mapSolution(VectorFiniteVolumeField &field)
{
    Teuchos::ArrayRCP<const Scalar> soln = x_->getData(0);
    for (const Cell &cell: field.grid().localActiveCells())
    {
        Vector2D &u = field(cell);
        u.x = soln[field.indexMap()->local(cell, 0)];
        u.y = soln[field.indexMap()->local(cell, 1)];
    }
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
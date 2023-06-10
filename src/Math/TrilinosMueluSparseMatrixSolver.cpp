#include "trilinos/MueLu_CreateTpetraPreconditioner.hpp"
#include "trilinos/MueLu_ParameterListInterpreter.hpp"
#include "trilinos/BelosSolverFactory.hpp"

#include "TrilinosMueluSparseMatrixSolver.h"

TrilinosMueluSparseMatrixSolver::TrilinosMueluSparseMatrixSolver(
    const Communicator &comm, const std::string &solverName)
    : TrilinosSparseMatrixSolver(comm) {
  typedef Belos::SolverFactory<Scalar, TpetraMultiVector, TpetraOperator>
      SolverFactory;

  belosParams_ = rcp(new Teuchos::ParameterList());
  mueluParams_ = rcp(new Teuchos::ParameterList());
  solver_ = SolverFactory().create(solverName, belosParams_);

  linearProblem_ = rcp(new LinearProblem());
  solver_->setProblem(linearProblem_);
}

void TrilinosMueluSparseMatrixSolver::setRank(int rank) {
  TrilinosSparseMatrixSolver::setRank(rank);
  linearProblem_->setOperator(mat_);
}

Scalar TrilinosMueluSparseMatrixSolver::solve() {
  precon_ = MueLu::CreateTpetraPreconditioner(
      Teuchos::rcp_static_cast<TpetraOperator>(mat_), *mueluParams_, coords_);

  linearProblem_->setProblem(x_, b_);
  linearProblem_->setLeftPrec(precon_);
  solver_->solve();

  return error();
}

void TrilinosMueluSparseMatrixSolver::setup(
    const boost::property_tree::ptree &parameters) {
  using namespace std;
  using namespace Teuchos;

  typedef Belos::SolverFactory<Scalar, TpetraMultiVector, TpetraOperator>
      SolverFactory;

  std::string belosParamFile = parameters.get<std::string>("belosParamFile");
  std::string mueluParamFile = parameters.get<std::string>("mueluParamFile");
  std::string solverName = parameters.get<std::string>("solver", "GMRES");

  belosParams_ = Teuchos::getParametersFromXmlFile("case/" + belosParamFile);
  mueluParams_ = Teuchos::getParametersFromXmlFile("case/" + mueluParamFile);
  solver_ = SolverFactory().create(solverName, belosParams_);

  linearProblem_ = rcp(new LinearProblem());
  solver_->setProblem(linearProblem_);
}

int TrilinosMueluSparseMatrixSolver::nIters() const {
  return solver_->getNumIters();
}

Scalar TrilinosMueluSparseMatrixSolver::error() const {
  return solver_->achievedTol();
}

void TrilinosMueluSparseMatrixSolver::printStatus(
    const std::string &msg) const {
  comm_.printf("%s %s iterations = %d, error = %lf.\n", msg.c_str(), "Krylov",
               nIters(), error());
}

void TrilinosMueluSparseMatrixSolver::setCoordinates(
    const std::vector<Point2D> &coordinates) {
  coords_ = Teuchos::rcp(new TpetraMultiVector(rangeMap_, 2));

  std::transform(coordinates.begin(), coordinates.end(),
                 coords_->getDataNonConst(0).begin(),
                 [](const Point2D &coord) { return coord.x; });

  std::transform(coordinates.begin(), coordinates.end(),
                 coords_->getDataNonConst(1).begin(),
                 [](const Point2D &coord) { return coord.y; });
}

void TrilinosMueluSparseMatrixSolver::setCoordinates(
    const std::vector<Point3D> &coordinates) {
  coords_ = Teuchos::rcp(new TpetraMultiVector(rangeMap_, 3));

  std::transform(coordinates.begin(), coordinates.end(),
                 coords_->getDataNonConst(0).begin(),
                 [](const Point3D &coord) { return coord.x; });

  std::transform(coordinates.begin(), coordinates.end(),
                 coords_->getDataNonConst(1).begin(),
                 [](const Point3D &coord) { return coord.y; });

  std::transform(coordinates.begin(), coordinates.end(),
                 coords_->getDataNonConst(2).begin(),
                 [](const Point3D &coord) { return coord.z; });
}

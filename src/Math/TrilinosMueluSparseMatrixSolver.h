#ifndef PHASE_TRILINOS_MUELU_SPARSE_MATRIX_SOLVER_H
#define PHASE_TRILINOS_MUELU_SPARSE_MATRIX_SOLVER_H

#include <BelosSolverManager.hpp>
#include <BelosTpetraAdapter.hpp>
#include <MueLu_TpetraOperator.hpp>

#include "Math/TrilinosSparseMatrixSolver.h"

#include "2D/Geometry/Point2D.h"
#include "3D/Geometry/Point3D.h"

class TrilinosMueluSparseMatrixSolver : public TrilinosSparseMatrixSolver {
public:
  TrilinosMueluSparseMatrixSolver(const Communicator &comm,
                                  const std::string &solverName = "TFQMR");

  Type type() const { return TRILINOS_MUELU; }

  void setRank(int rank);

  Scalar solve();

  void setup(const boost::property_tree::ptree &parameters);

  int nIters() const;

  Scalar error() const;

  void printStatus(const std::string &msg) const;

  void setCoordinates(const std::vector<Point2D> &coordinates);

  void setCoordinates(const std::vector<Point3D> &coordinates);

private:
  typedef Belos::LinearProblem<Scalar, TpetraMultiVector, TpetraOperator>
      LinearProblem;
  typedef Belos::SolverManager<Scalar, TpetraMultiVector, TpetraOperator>
      Solver;
  typedef MueLu::TpetraOperator<Scalar, Index, Index> Preconditioner;

  Teuchos::RCP<Teuchos::ParameterList> belosParams_, mueluParams_;

  Teuchos::RCP<TpetraMultiVector> coords_;

  Teuchos::RCP<LinearProblem> linearProblem_;

  Teuchos::RCP<Solver> solver_;

  Teuchos::RCP<Preconditioner> precon_;
};

#endif

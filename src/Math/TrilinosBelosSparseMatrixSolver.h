#ifndef TRILINOS_BELOS_SPARSE_MATRIX_SOLVER_H
#define TRILINOS_BELOS_SPARSE_MATRIX_SOLVER_H

#include "BelosSolverManager.hpp"
#include "BelosTpetraAdapter.hpp"
#include "Ifpack2_Preconditioner.hpp"

#include "TrilinosSparseMatrixSolver.h"

class TrilinosBelosSparseMatrixSolver : public TrilinosSparseMatrixSolver {
public:
  TrilinosBelosSparseMatrixSolver(const Communicator &comm);

  Type type() const { return TRILINOS_BELOS; }

  void setRank(int rank);

  Scalar solve();

  void setup(const boost::property_tree::ptree &parameters);

  int nIters() const;

  Scalar error() const;

  void printStatus(const std::string &msg) const;

private:
  typedef Belos::LinearProblem<Scalar, TpetraMultiVector, TpetraOperator>
      LinearProblem;
  typedef Belos::SolverManager<Scalar, TpetraMultiVector, TpetraOperator>
      Solver;
  typedef Ifpack2::Preconditioner<Scalar, Index, Index> Preconditioner;

  //- Types
  std::string precType_;

  //- Parameters
  Teuchos::RCP<Teuchos::ParameterList> belosParams_, ifpackParams_;

  //- Solver data structures
  Teuchos::RCP<LinearProblem> linearProblem_;
  Teuchos::RCP<Solver> solver_;
  Teuchos::RCP<Preconditioner> precon_;
};

#endif

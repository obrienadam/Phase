#ifndef PHASE_SPARSE_MATRIX_SOLVER_FACTORY_H
#define PHASE_SPARSE_MATRIX_SOLVER_FACTORY_H

#include "System/Communicator.h"

#include "SparseMatrixSolver.h"

class SparseMatrixSolverFactory {
public:
  enum Type { EIGEN, TRILINOS_BELOS, TRILINOS_AMESOS2, TRILINOS_MUELU };

  std::shared_ptr<SparseMatrixSolver> create(Type type,
                                             const Communicator &comm) const;

  std::shared_ptr<SparseMatrixSolver> create(const std::string &type,
                                             const Communicator &comm) const;

protected:
};

#endif

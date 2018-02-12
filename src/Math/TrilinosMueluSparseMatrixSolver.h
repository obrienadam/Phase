#ifndef TRILINOS_MUELU_SPARSE_MATRIX_SOLVER_H
#define TRILINOS_MUELU_SPARSE_MATRIX_SOLVER_H

#include <MueLu_TpetraOperator.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosSolverManager.hpp>

#include "TrilinosSparseMatrixSolver.h"

class TrilinosMueluSparseMatrixSolver : public TrilinosSparseMatrixSolver
{
public:

    TrilinosMueluSparseMatrixSolver(const Communicator &comm,
                                    const std::weak_ptr<const FiniteVolumeGrid2D> &grid,
                                    const std::string &solverName = "BiCGSTAB");

    void setRank(int rank);

    Scalar solve();

    void setup(const boost::property_tree::ptree &parameters);

    int nIters() const;

    Scalar error() const;

    void printStatus(const std::string &msg) const;

private:

    typedef Belos::LinearProblem<Scalar, TpetraMultiVector, TpetraOperator> LinearProblem;
    typedef Belos::SolverManager<Scalar, TpetraMultiVector, TpetraOperator> Solver;
    typedef MueLu::TpetraOperator<Scalar, Index, Index> Preconditioner;

    Teuchos::RCP<Teuchos::ParameterList> belosParams_, mueluParams_;

    Teuchos::RCP<TpetraMultiVector> coords_;

    Teuchos::RCP<LinearProblem> linearProblem_;

    Teuchos::RCP<Solver> solver_;

    Teuchos::RCP<Preconditioner> precon_;

    std::weak_ptr<const FiniteVolumeGrid2D> grid_;
};

#endif
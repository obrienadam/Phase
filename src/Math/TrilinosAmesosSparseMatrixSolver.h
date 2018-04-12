#ifndef TRILINOS_AMESOS_SPARSE_MATRIX_SOLVER_H
#define TRILINOS_AMESOS_SPARSE_MATRIX_SOLVER_H

#include <Amesos2.hpp>

#include "TrilinosSparseMatrixSolver.h"

class TrilinosAmesosSparseMatrixSolver: public TrilinosSparseMatrixSolver
{
public:

    TrilinosAmesosSparseMatrixSolver(const Communicator &comm, const std::string &solverName = "klu2");

    Type type() const
    { return TRILINOS_AMESOS2; }

    Scalar solve();

    void setup(const boost::property_tree::ptree& parameters);

    int nIters() const
    { return 1; }

    Scalar error() const
    { return 0.; }

private:

    typedef Amesos2::Solver<TpetraCrsMatrix, TpetraMultiVector> Solver;

    std::string solverName_;

    Teuchos::RCP<Teuchos::ParameterList> amesos2Params_;

    Teuchos::RCP<Solver> solver_;
};

#endif
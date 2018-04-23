#include "System/Exception.h"

#include "SparseMatrixSolverFactory.h"

#include "EigenSparseMatrixSolver.h"
#include "TrilinosBelosSparseMatrixSolver.h"
#include "TrilinosAmesosSparseMatrixSolver.h"
#include "TrilinosMueluSparseMatrixSolver.h"

std::shared_ptr<SparseMatrixSolver> SparseMatrixSolverFactory::create(Type type, const Communicator &comm) const
{
    switch (type)
    {
        case EIGEN:
            return std::make_shared<EigenSparseMatrixSolver>();
        case TRILINOS_BELOS:
            return std::make_shared<TrilinosBelosSparseMatrixSolver>(comm);
        case TRILINOS_AMESOS2:
            return std::make_shared<TrilinosAmesosSparseMatrixSolver>(comm);
        case TRILINOS_MUELU:
            return std::make_shared<TrilinosMueluSparseMatrixSolver>(comm);
        default:
            return nullptr;
    }
}

std::shared_ptr<SparseMatrixSolver> SparseMatrixSolverFactory::create(const std::string &type,
                                                                      const Communicator &comm) const
{
    if (type == "eigen" || type == "eigen3")
        return create(EIGEN, comm);
    else if (type == "belos")
        return create(TRILINOS_BELOS, comm);
    else if (type == "amesos" || type == "amesos2")
        return create(TRILINOS_AMESOS2, comm);
    else if (type == "muelu")
        return create(TRILINOS_MUELU, comm);
    else
        throw Exception("SparseMatrixSolverFactory", "create", "bad solver type \"" + type + "\".");
}
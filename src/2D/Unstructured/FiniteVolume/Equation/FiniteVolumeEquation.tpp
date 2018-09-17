#include <stdio.h>

#include "System/Exception.h"

#include "FiniteVolumeEquation.h"

#include "Math/SparseMatrixSolverFactory.h"
#include "Math/TrilinosMueluSparseMatrixSolver.h"

template<class T>
FiniteVolumeEquation<T>::FiniteVolumeEquation(const Input &input,
                                              FiniteVolumeField<T> &field,
                                              const std::string &name, int nnz)
    :
      FiniteVolumeEquation<T>::FiniteVolumeEquation(field, name, nnz)
{
    configureSparseSolver(input, field.grid()->comm());
}

template<class T>
FiniteVolumeEquation<T> &FiniteVolumeEquation<T>::operator =(const FiniteVolumeEquation<T> &rhs)
{
    return operator =(static_cast<const CrsEquation&>(rhs));
}

template<class T>
FiniteVolumeEquation<T> &FiniteVolumeEquation<T>::operator =(FiniteVolumeEquation<T> &&rhs)
{
    return operator =(static_cast<CrsEquation&&>(rhs));
}

template<class T>
FiniteVolumeEquation<T> &FiniteVolumeEquation<T>::operator =(const CrsEquation &rhs)
{
    if(this != &rhs)
        CrsEquation::operator =(rhs);
    return *this;
}

template<class T>
FiniteVolumeEquation<T> &FiniteVolumeEquation<T>::operator =(CrsEquation &&rhs)
{
    CrsEquation::operator =(rhs);
    return *this;
}

template<class T>
void FiniteVolumeEquation<T>::configureSparseSolver(const Input &input, const Communicator &comm)
{
    std::string lib = input.caseInput().get<std::string>("LinearAlgebra." + name + ".lib");
    boost::algorithm::to_lower(lib);

    solver_ = SparseMatrixSolverFactory().create(lib, comm);

    if (comm.nProcs() > 1 && !solver_->supportsMPI())
        throw Exception("FiniteVolumeEquation<T>", "configureSparseSolver", "equation \"" + name + "\", lib \"" + lib +
                        "\" does not support multiple processes in its current configuration.");

    solver_->setup(input.caseInput().get_child("LinearAlgebra." + name));

    comm.printf("Initialized sparse matrix solver for equation \"%s\" using lib%s.\n", name.c_str(), lib.c_str());
}

template<class T>
Scalar FiniteVolumeEquation<T>::solve()
{
    if (!solver_)
        throw Exception("FiniteVolumeEquation<T>", "solve",
                        "must allocate a SparseMatrixSolver object before attempting to solve.");

    solver_->setRank(getRank());
    solver_->set(rowPtr_, colInd_, vals_);
    solver_->setRhs(-rhs_);

    if (solver_->type() == SparseMatrixSolver::TRILINOS_MUELU)
        std::static_pointer_cast<TrilinosMueluSparseMatrixSolver>(solver_)->setCoordinates(
                    field_.grid()->localCells().coordinates());

    solver_->solve();

    mapFromSparseSolver();

    solver_->printStatus("FiniteVolumeEquation " + name + ":");

    return solver_->error();
}

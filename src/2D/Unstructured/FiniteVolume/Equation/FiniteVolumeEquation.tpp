#include <stdio.h>

#include "System/Exception.h"

#include "FiniteVolumeEquation.h"

#include "Math/SparseMatrixSolverFactory.h"
#include "Math/TrilinosMueluSparseMatrixSolver.h"

template<class T>
FiniteVolumeEquation<T>::FiniteVolumeEquation(const Input &input,
                                              FiniteVolumeField<T> &field,
                                              const std::string &name)
    :
      FiniteVolumeEquation<T>::FiniteVolumeEquation(field, name)
{
    configureSparseSolver(input, field.grid()->comm());
}

template<class T>
FiniteVolumeEquation<T> &FiniteVolumeEquation<T>::operator =(const FiniteVolumeEquation<T> &rhs)
{
    return operator =(static_cast<const Equation&>(rhs));
}

template<class T>
FiniteVolumeEquation<T> &FiniteVolumeEquation<T>::operator =(FiniteVolumeEquation<T> &&rhs)
{
    return operator =(static_cast<Equation&&>(rhs));
}

template<class T>
FiniteVolumeEquation<T> &FiniteVolumeEquation<T>::operator =(const Equation &rhs)
{
    if(this != &rhs)
        Equation::operator =(rhs);
    return *this;
}

template<class T>
FiniteVolumeEquation<T> &FiniteVolumeEquation<T>::operator =(Equation &&rhs)
{
    Equation::operator =(rhs);
    return *this;
}

template<class T>
void FiniteVolumeEquation<T>::configureSparseSolver(const Input &input, const Communicator &comm)
{
    std::string lib = input.caseInput().get<std::string>("LinearAlgebra." + name + ".lib");
    boost::algorithm::to_lower(lib);

    _spSolver = SparseMatrixSolverFactory().create(lib, comm);

    if (comm.nProcs() > 1 && !_spSolver->supportsMPI())
        throw Exception("FiniteVolumeEquation<T>", "configureSparseSolver", "equation \"" + name + "\", lib \"" + lib +
                        "\" does not support multiple processes in its current configuration.");

    _spSolver->setup(input.caseInput().get_child("LinearAlgebra." + name));

    comm.printf("Initialized sparse matrix solver for equation \"%s\" using lib%s.\n", name.c_str(), lib.c_str());
}

template<class T>
Scalar FiniteVolumeEquation<T>::solve()
{
    if (!_spSolver)
        throw Exception("FiniteVolumeEquation<T>", "solve",
                        "must allocate a SparseMatrixSolver object before attempting to solve.");

    _spSolver->setRank(getRank());
    _spSolver->set(_coeffs);
    _spSolver->setRhs(-_rhs);

    if (_spSolver->type() == SparseMatrixSolver::TRILINOS_MUELU)
        std::static_pointer_cast<TrilinosMueluSparseMatrixSolver>(_spSolver)->setCoordinates(
                    field_.grid()->localCells().coordinates());

    _spSolver->solve();

    mapFromSparseSolver();

    _spSolver->printStatus("FiniteVolumeEquation " + name + ":");

    return _spSolver->error();
}

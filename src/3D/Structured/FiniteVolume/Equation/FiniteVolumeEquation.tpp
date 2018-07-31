#include "System/Exception.h"

#include "FiniteVolumeEquation.h"
#include "Math/SparseMatrixSolverFactory.h"

template<class T>
FiniteVolumeEquation<T>::FiniteVolumeEquation(const std::string &name, Field<T> &field)
    :
      FiniteVolumeEquation(field)
{
    _name = name;
}

template<class T>
FiniteVolumeEquation<T>::FiniteVolumeEquation(const std::string &name, const Input &input, Field<T> &field)
    :
      FiniteVolumeEquation(name, field)
{
    configureSparseSolver(input);
}

template<class T>
FiniteVolumeEquation<T>& FiniteVolumeEquation<T>::operator=(const Equation &rhs)
{
    _coeffs = rhs.coeffs();
    _rhs = rhs.rhs();
    return *this;
}

template<class T>
void FiniteVolumeEquation<T>::configureSparseSolver(const Input &input)
{
    std::string lib = input.caseInput().get<std::string>("LinearAlgebra." + _name + ".lib");
    boost::algorithm::to_lower(lib);

    _spSolver = SparseMatrixSolverFactory().create(lib, _field.grid()->comm());

    if (_field.grid()->comm().nProcs() > 1 && !_spSolver->supportsMPI())
        throw Exception("FiniteVolumeEquation<T>",
                        "configureSparseSolver",
                        "equation \"" + _name + "\", lib \""
                        + lib + "\" does not support multiple processes in its current configuration.");

    _spSolver->setup(input.caseInput().get_child("LinearAlgebra." + _name));

    _field.grid()->comm().printf("Initialized sparse matrix solver for equation \"%s\" using lib%s.\n",
                                 _name.c_str(),
                                 lib.c_str());
}

template<class T>
Scalar FiniteVolumeEquation<T>::solve()
{
    Equation::solve();
    for (const Cell &cell: _field.grid()->cells())
        _field(cell) = _spSolver->x(cell.id());
}

#include "Solver.h"

Solver::Solver(const Input &input, const std::shared_ptr<StructuredGrid2D> &grid)
    :
      _grid(grid)
{

}

void Solver::setInitialConditions(const Input &input)
{

}

void Solver::setInitialConditions(const CommandLine &cl, const Input &input)
{

}

template<>
ScalarField& Solver::addField(const std::string &name)
{
    auto insert = _scalarFields.insert(std::make_shared<ScalarField>(name, _grid, true, true));
    return *(*insert.first);
}

template<>
VectorField& Solver::addField(const std::string &name)
{
    auto insert = _vectorFields.insert(std::make_shared<VectorField>(name, _grid, true, true));
    return *(*insert.first);
}

template<>
Field<Scalar> &Solver::addField(const std::string &name, const Input &input)
{
    auto insert = _scalarFields.insert(std::make_shared<ScalarField>(name, _grid, input, true, true));
    return *(*insert.first);
}

template<>
Field<Vector2D> &Solver::addField(const std::string &name, const Input &input)
{
    auto insert = _vectorFields.insert(std::make_shared<VectorField>(name, _grid, input, true, true));
    return *(*insert.first);
}

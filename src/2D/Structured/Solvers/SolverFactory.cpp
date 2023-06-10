#include "System/Exception.h"

#include "Poisson.h"
#include "SolverFactory.h"

std::shared_ptr<Solver>
SolverFactory::create(Type type, const Input &input,
                      const std::shared_ptr<StructuredGrid2D> &grid) {
  switch (type) {
  case POISSON:
    return std::make_shared<Poisson>(input, grid);
  }
}

std::shared_ptr<Solver>
SolverFactory::create(std::string type, const Input &input,
                      const std::shared_ptr<StructuredGrid2D> &grid) {
  std::transform(type.begin(), type.end(), type.begin(), ::tolower);

  if (type == "poisson")
    return create(POISSON, input, grid);
  else
    throw Exception("SolverFactory", "create",
                    "solver type \"" + type + "\" not recognized.");
}

std::shared_ptr<Solver>
SolverFactory::create(const Input &input,
                      const std::shared_ptr<StructuredGrid2D> &grid) {
  return create(input.caseInput().get<std::string>("Solver.type"), input, grid);
}

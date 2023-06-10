#include "SolverFactory.h"
#include "Poisson.h"

std::shared_ptr<Solver>
SolverFactory::create(SolverFactory::SolverType type, const Input &input,
                      const std::shared_ptr<const StructuredGrid3D> &grid) {
  switch (type) {
  case POISSON:
    return std::make_shared<Poisson>(input, grid);
  case FRACTIONAL_STEP:
    return nullptr;
  }
}

std::shared_ptr<Solver>
SolverFactory::create(std::string type, const Input &input,
                      const std::shared_ptr<const StructuredGrid3D> &grid) {
  std::transform(type.begin(), type.end(), type.begin(),
                 [](unsigned char c) { return std::tolower(c); });

  if (type == "poisson")
    return create(POISSON, input, grid);
  else if (type == "fractional step")
    return create(FRACTIONAL_STEP, input, grid);
}

std::shared_ptr<Solver>
SolverFactory::create(const Input &input,
                      const std::shared_ptr<const StructuredGrid3D> &grid) {
  return create(input.caseInput().get<std::string>("Solver.type"), grid);
}

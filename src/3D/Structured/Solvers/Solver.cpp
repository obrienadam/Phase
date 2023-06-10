#include "Solver.h"

Solver::Solver(const Input &input,
               const std::shared_ptr<const StructuredGrid3D> &grid)
    : _grid(grid) {}

Scalar Solver::getStartTime() const {}

std::string Solver::info() const { return "Unknown structured 3D solver"; }

void Solver::setInitialConditions(const Input &input) {}

void Solver::setInitialConditions(const CommandLine &cl, const Input &input) {}

int Solver::printf(const char *format, ...) const {
  int n = 0;

  if (_grid->comm().isMainProc()) {
    va_list argsPtr;
    va_start(argsPtr, format);
    n = vfprintf(stdout, format, argsPtr);
    va_end(argsPtr);
  }

  return n;
}

template <>
const std::shared_ptr<ScalarField> &Solver::addField(const std::string &name,
                                                     const Input &input) {
  return _scalarFields[name] =
             std::make_shared<ScalarField>(name, input, _grid, true, true);
}

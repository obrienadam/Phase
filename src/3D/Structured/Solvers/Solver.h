#ifndef PHASE_SOLVER_H
#define PHASE_SOLVER_H

#include "Structured/FiniteVolume/Field/ScalarField.h"
#include "Structured/StructuredGrid3D/StructuredGrid3D.h"
#include "System/SolverInterface.h"

class Solver : public SolverInterface {
public:
  Solver(const Input &input,
         const std::shared_ptr<const StructuredGrid3D> &grid);

  virtual Scalar getStartTime() const override;

  virtual std::string info() const override;

  virtual void initialize() override {}

  void setInitialConditions(const Input &input) override;

  void setInitialConditions(const CommandLine &cl, const Input &input) override;

  virtual Scalar maxTimeStep() const override {
    return std::numeric_limits<Scalar>::infinity();
  }

  virtual Scalar computeMaxTimeStep(Scalar maxCo,
                                    Scalar timeStep) const override {
    return std::numeric_limits<Scalar>::infinity();
  }

  virtual Scalar solve(Scalar timeStep) = 0;

  int printf(const char *format, ...) const override;

  const std::shared_ptr<const StructuredGrid3D> &grid() const { return _grid; }

  const Communicator &comm() const override { return _grid->comm(); }

  //- Add fields

  template <class T>
  const std::shared_ptr<Field<T>> &addField(const std::string &name,
                                            const Input &input);

  //- Get fields

  const std::unordered_map<std::string, std::shared_ptr<ScalarField>> &
  scalarFields() const {
    return _scalarFields;
  }

protected:
  std::shared_ptr<const StructuredGrid3D> _grid;

  std::unordered_map<std::string, std::shared_ptr<ScalarField>> _scalarFields;
};

#endif

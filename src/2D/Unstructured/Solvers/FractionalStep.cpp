#include "FiniteVolume/Discretization/Divergence.h"
#include "FiniteVolume/Discretization/Laplacian.h"
#include "FiniteVolume/Discretization/Source.h"
#include "FiniteVolume/Discretization/TimeDerivative.h"

#include "FractionalStep.h"

FractionalStep::FractionalStep(
    const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
    : Solver(input, grid), fluid_(std::make_shared<CellGroup>("fluid")),
      u_(*addField<Vector2D>(input, "u", fluid_)),
      p_(*addField<Scalar>(input, "p", fluid_)),
      co_(*addField<Scalar>(input, "co", fluid_)),
      gradP_(*std::static_pointer_cast<ScalarGradient>(
          addField<Vector2D>(std::make_shared<ScalarGradient>(p_, fluid_)))),
      gradU_(*std::static_pointer_cast<JacobianField>(
          addField<Tensor2D>(std::make_shared<JacobianField>(u_, fluid_)))),
      uEqn_(input, u_, "uEqn"), pEqn_(input, p_, "pEqn") {
  fluid_->add(grid_->localCells());
  rho_ = input.caseInput().get<Scalar>("Properties.rho", 1);
  mu_ = input.caseInput().get<Scalar>("Properties.mu", 1);
  g_ = input.caseInput().get<std::string>("Properties.g", "(0,0)");
}

void FractionalStep::initialize() {
  u_.interpolateFaces();
  p_.setBoundaryFaces();
}

std::string FractionalStep::info() const {
  return "Fractional-step\n"
         "A simple 1-step fractional-step projection method\n"
         "May not produce accurate results near boundaries\n";
}

Scalar FractionalStep::solve(Scalar timeStep) {
  solveUEqn(timeStep);
  solvePEqn(timeStep);
  correctVelocity(timeStep);

  grid_->comm().printf("Max divergence error = %.4e\n",
                       grid_->comm().max(maxDivergenceError()));
  grid_->comm().printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

  return 0;
}

Scalar FractionalStep::maxCourantNumber(Scalar timeStep) const {
  Scalar maxCo = 0;

  for (const Cell &cell : *fluid_) {
    Scalar co = 0.;

    for (const InteriorLink &nb : cell.neighbours())
      co += std::max(dot(u_(nb.face()), nb.outwardNorm()), 0.);

    for (const BoundaryLink &bd : cell.boundaries())
      co += std::max(dot(u_(bd.face()), bd.outwardNorm()), 0.);

    co *= timeStep / cell.volume();
    co_(cell) = co;
    maxCo = std::max(co, maxCo);
  }

  return grid_->comm().max(maxCo);
}

Scalar FractionalStep::computeMaxTimeStep(Scalar maxCo,
                                          Scalar prevTimeStep) const {
  Scalar co = maxCourantNumber(prevTimeStep);
  Scalar lambda1 = 0.1, lambda2 = 1.2;

  return grid_->comm().min(
      std::min(std::min(maxCo / co * prevTimeStep,
                        (1 + lambda1 * maxCo / co) * prevTimeStep),
               std::min(lambda2 * prevTimeStep, maxTimeStep_)));
}

Scalar FractionalStep::solveUEqn(Scalar timeStep) {
  u_.savePreviousTimeStep(timeStep, 1);

  uEqn_ = (fv::ddt(u_, timeStep) + fv::div(u_, u_, 0.) ==
           fv::laplacian(mu_ / rho_, u_, 0.5) - src::src(gradP_));

  Scalar error = uEqn_.solve();

  for (const Cell &cell : *fluid_)
    u_(cell) += timeStep * gradP_(cell);

  grid_->sendMessages(u_);
  u_.interpolateFaces();

  return error;
}

Scalar FractionalStep::solvePEqn(Scalar timeStep) {
  pEqn_ = (fv::laplacian(timeStep, p_) == src::div(u_));

  Scalar error = pEqn_.solve();
  grid_->sendMessages(p_);
  p_.setBoundaryFaces();

  //- Gradient
  gradP_.compute(*fluid_);

  return error;
}

void FractionalStep::correctVelocity(Scalar timeStep) {
  for (const Cell &cell : *fluid_)
    u_(cell) -= timeStep * gradP_(cell);

  grid_->sendMessages(u_); //- Necessary

  for (const Face &face : grid_->faces())
    u_(face) -= timeStep * gradP_(face);
}

Scalar FractionalStep::maxDivergenceError() {
  Scalar maxError = 0.;

  for (const Cell &cell : *fluid_) {
    Scalar div = 0.;

    for (const InteriorLink &nb : cell.neighbours())
      div += dot(u_(nb.face()), nb.outwardNorm());

    for (const BoundaryLink &bd : cell.boundaries())
      div += dot(u_(bd.face()), bd.outwardNorm());

    maxError = fabs(div) > maxError ? div : maxError;
  }

  return grid_->comm().max(maxError);
}

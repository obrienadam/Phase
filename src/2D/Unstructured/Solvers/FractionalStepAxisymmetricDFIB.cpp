#include "FiniteVolume/Discretization/AxisymmetricDivergence.h"
#include "FiniteVolume/Discretization/AxisymmetricExplicitDivergence.h"
#include "FiniteVolume/Discretization/AxisymmetricLaplacian.h"
#include "FiniteVolume/Discretization/AxisymmetricSource.h"
#include "FiniteVolume/Discretization/AxisymmetricStressTensor.h"
#include "FiniteVolume/Discretization/AxisymmetricTimeDerivative.h"
#include "FiniteVolume/ImmersedBoundary/DirectForcingImmersedBoundary.h"
#include "FiniteVolume/Motion/SolidBodyMotion.h"

#include "FractionalStepAxisymmetricDFIB.h"

FractionalStepAxisymmetricDFIB::FractionalStepAxisymmetricDFIB(
    const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
    : FractionalStepAxisymmetric(input, grid),
      fib_(*addField<Vector2D>("fb", fluid_)), fibEqn_(input, fib_, "fbEqn") {
  ib_ = std::make_shared<DirectForcingImmersedBoundary>(input, grid, fluid_);
  addField<int>(ib_->cellStatus());

  for (auto &ibObj : *ib_) {
    auto motion = std::dynamic_pointer_cast<SolidBodyMotion>(ibObj->motion());

    if (motion)
      motion->setMotionConstraint(Vector2D(0., 1.));
  }
}

Scalar FractionalStepAxisymmetricDFIB::solve(Scalar timeStep) {
  grid_->comm().printf("Updating IB positions...\n");
  ib_->updateIbPositions(timeStep);
  ib_->updateCells();

  solveUEqn(timeStep);
  solvePEqn(timeStep);
  correctVelocity(timeStep);
  computeIbForces(timeStep);

  grid_->comm().printf("Max divergence error = %.4e\n",
                       grid_->comm().max(maxDivergenceError()));
  grid_->comm().printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

  return 0.;
}

std::shared_ptr<const ImmersedBoundary>
FractionalStepAxisymmetricDFIB::ib() const {
  return ib_;
}

Scalar FractionalStepAxisymmetricDFIB::solveUEqn(Scalar timeStep) {
  u_.savePreviousTimeStep(timeStep, 2);
  uEqn_ = (axi::ddt(u_, timeStep) + axi::dive(u_, u_, 0.5) ==
           axi::laplacian(mu_ / rho_, u_, 0.5) - axi::src::src(gradP_));

  Scalar error = uEqn_.solve();
  u_.sendMessages();

  fibEqn_ = ib_->computeForcingTerm(u_, timeStep, fib_);
  fibEqn_.solve();
  fib_.sendMessages();

  for (const Cell &c : *fluid_) {
    u_(c) += timeStep * fib_(c);
    u_(c) += timeStep * gradP_(c);
  }

  u_.sendMessages();
  u_.interpolateFaces();

  return error;
}

void FractionalStepAxisymmetricDFIB::computeIbForces(Scalar timeStep) {
  for (auto &ibObj : *ib_) {
    if (ibObj->shape().type() != Shape2D::CIRCLE ||
        ibObj->shape().centroid().x != 0.)
      throw Exception("FractionalStepAxisymmetricDFIB", "computeIbForces",
                      "oncly circles centered at r = 0 supported.");

    Vector2D fh(0., 0.);

    for (const Cell &c : ibObj->cells()) {
      fh += (u_(c) - u_.oldField(0)(c)) * c.polarVolume() / timeStep;

      for (const InteriorLink &nb : c.neighbours()) {
        Scalar flux0 =
            dot(u_.oldField(0)(nb.face()), nb.polarOutwardNorm()) / 2.;
        Scalar flux1 =
            dot(u_.oldField(1)(nb.face()), nb.polarOutwardNorm()) / 2.;
        fh += std::max(flux0, 0.) * u_.oldField(0)(c) +
              std::min(flux0, 0.) * u_.oldField(0)(nb.cell()) +
              std::max(flux1, 0.) * u_.oldField(1)(c) +
              std::min(flux1, 0.) * u_.oldField(1)(nb.cell());
      }

      for (const BoundaryLink &bd : c.boundaries()) {
        Scalar flux0 =
            dot(u_.oldField(0)(bd.face()), bd.polarOutwardNorm()) / 2.;
        Scalar flux1 =
            dot(u_.oldField(1)(bd.face()), bd.polarOutwardNorm()) / 2.;
        fh += std::max(flux0, 0.) * u_.oldField(0)(c) +
              std::min(flux0, 0.) * u_.oldField(0)(bd.face()) +
              std::max(flux1, 0.) * u_.oldField(1)(c) +
              std::min(flux1, 0.) * u_.oldField(1)(bd.face());
      }

      fh -= fib_(c) * c.polarVolume();
    }

    fh = grid_->comm().sum(2. * M_PI * rho_ * fh);

    //- Assume spherical
    const Circle &circ = static_cast<const Circle &>(ibObj->shape());
    Scalar vol = 4. / 3. * M_PI * std::pow(circ.radius(), 3);
    Vector2D fw = (ibObj->rho - rho_) * vol * g_;

    // std::cout << "Cd = " << 2. * fh.y / (rho_ * M_PI *
    // std::pow(circ.radius(), 2)) << "\n";

    ibObj->applyForce((fh + fw) * ibObj->mass() / (ibObj->rho * vol));
  }
}

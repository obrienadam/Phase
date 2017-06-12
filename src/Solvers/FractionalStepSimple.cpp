#include "FractionalStepSimple.h"
#include "FaceInterpolation.h"
#include "GradientEvaluation.h"
#include "Source.h"
#include "SeoMittal.h"

FractionalStepSimple::FractionalStepSimple(const Input &input,
                                           const Communicator &comm,
                                           std::shared_ptr<FiniteVolumeGrid2D> &grid)
    :
      Solver(input, comm, grid),
      u(addVectorField(input, "u")),
      gradP(addVectorField("gradP")),
      p(addScalarField(input, "p")),
      uEqn_(input, comm, u, "uEqn"),
      pEqn_(input, comm, p, "pEqn"),
      fluid_(grid->createCellZone("fluid"))
{
    rho_ = input.caseInput().get<Scalar>("Properties.rho", 1);
    mu_ = input.caseInput().get<Scalar>("Properties.mu", 1);

    //- All active cells to fluid cells
    fluid_.add(grid_->localActiveCells());

    //- Create ib zones if any. Will also update local/global indices
    ib_.initCellZones(fluid_);
}

void FractionalStepSimple::initialize()
{

}

std::string FractionalStepSimple::info() const
{
    return Solver::info()
            + "Fractional-step (Simple)\n"
            + "A simple 1-step fractional-step projection method\n"
            + "May not produce accurate results near boundaries\n";
}

Scalar FractionalStepSimple::solve(Scalar timeStep)
{
    solveUEqn(timeStep);

    ib_.clearFreshCells();

    solvePEqn(timeStep);
    correctVelocity(timeStep);

    ib_.computeForce(ScalarFiniteVolumeField(grid(), "rho", rho_),
                     ScalarFiniteVolumeField(grid(), "mu", mu_), u, p);

    ib_.update(timeStep);

    comm_.printf("Max divergence error = %.4e\n", comm_.max(seo::maxDivergence(ib_, u)));
    comm_.printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

    return 0;
}

Scalar FractionalStepSimple::maxCourantNumber(Scalar timeStep) const
{
    Scalar maxCo = 0;

    for (const Face &face: grid_->interiorFaces())
    {
        Vector2D sf = face.outwardNorm(face.lCell().centroid());
        Vector2D rc = face.rCell().centroid() - face.lCell().centroid();

        maxCo = std::max(maxCo, fabs(dot(u(face), sf) / dot(rc, sf)));
    }

    return comm_.max(maxCo * timeStep);
}

Scalar FractionalStepSimple::computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
{
    Scalar co = maxCourantNumber(prevTimeStep);
    Scalar lambda1 = 0.1, lambda2 = 1.2;

    return comm_.min(
                std::min(
                    std::min(maxCo / co * prevTimeStep, (1 + lambda1 * maxCo / co) * prevTimeStep),
                    std::min(lambda2 * prevTimeStep, maxTimeStep_)
                    ));
}

Scalar FractionalStepSimple::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    uEqn_ = (fv::ddt(u, timeStep) + fv::div(u, u) + ib_.bcs(u) == fv::laplacian(mu_/rho_, u));
    Scalar error = uEqn_.solve();

    u.interpolateFaces([](const Face& f){
        Scalar l1 = (f.lCell().centroid() - f.centroid()).mag();
        Scalar l2 = (f.rCell().centroid() - f.centroid()).mag();

        return l2/(l1 + l2);
    });

    return error;
}

Scalar FractionalStepSimple::solvePEqn(Scalar timeStep)
{
    pEqn_ = (seo::laplacian(ib_, rho_, timeStep, p) == seo::div(ib_, u));
    //pEqn_ = (fv::laplacian(timeStep/rho_, p) + ib_.bcs(p) == source::div(u));
    Scalar error = pEqn_.solve();

    p.interpolateFaces([](const Face& f){
        Scalar l1 = (f.lCell().centroid() - f.centroid()).mag();
        Scalar l2 = (f.rCell().centroid() - f.centroid()).mag();

        return l2/(l1 + l2);
    });

    return error;
}

void FractionalStepSimple::correctVelocity(Scalar timeStep)
{
    seo::correct(ib_, rho_, p, gradP, u, timeStep);
    return;

    for(const Cell& cell: fluid_)
        u(cell) -= timeStep/rho_*gradP(cell);

    for(const Face& face: grid_->interiorFaces())
        u(face) -= timeStep/rho_*gradP(face);

    u.setBoundaryFaces();
}

Scalar FractionalStepSimple::maxDivergenceError()
{
    Scalar maxError = 0.;

    for (const Cell &cell: fluid_)
    {
        Scalar div = 0.;

        for (const InteriorLink &nb: cell.neighbours())
            div += dot(u(nb.face()), nb.outwardNorm());

        for (const BoundaryLink &bd: cell.boundaries())
            div += dot(u(bd.face()), bd.outwardNorm());

        if (fabs(div) > maxError)
            maxError = fabs(div);
    }

    return comm_.max(maxError);
}



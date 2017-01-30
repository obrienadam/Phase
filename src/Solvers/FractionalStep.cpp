#include "FractionalStep.h"
#include "CrankNicolson.h"
#include "GradientEvaluation.h"
#include "FaceInterpolation.h"
#include "SourceEvaluation.h"

FractionalStep::FractionalStep(const Input &input, const Communicator &comm, FiniteVolumeGrid2D& grid)
    :
      Solver(input, comm, grid),
      u(addVectorField(input, "u")),
      gradP(addVectorField("gradP")),
      gradDp(addVectorField("gradDp")),
      p(addScalarField(input, "p")),
      dp(addScalarField("dp")),
      rho(addScalarField("rho")),
      mu(addScalarField("mu")),
      divUStar(addScalarField("uStar")),
      uEqn_(input, comm, u, "uEqn"),
      pEqn_(input, comm, dp, "pEqn")
{
    rho.fill(input.caseInput().get<Scalar>("Properties.rho", 1.));
    mu.fill(input.caseInput().get<Scalar>("Properties.mu", 1.));

    rho.savePreviousTimeStep(0., 1);
    mu.savePreviousTimeStep(0., 1);

    dp.copyBoundaryTypes(p);

    //- All active cells to fluid cells
    grid_.createCellZone("fluid", grid_.getCellIds(grid_.activeCells()));
    ib_.initCellZones();
    grid_.computeOrdering(comm_);

    volumeIntegrators_ = VolumeIntegrator::initVolumeIntegrators(input, *this);
    forceIntegrators_ = ForceIntegrator::initForceIntegrators(input, p, rho, mu, u);
}

std::string FractionalStep::info() const
{
    return Solver::info()
            + "Type: 2nd order fractional-step\n"
            + "Advection time-marching: Crank-Nicolson\n"
            + "Diffusion time-marching: Adams-Bashforth\n";
}

Scalar FractionalStep::solve(Scalar timeStep)
{
    solveUEqn(timeStep);
    solvePEqn(timeStep); comm_.barrier();
    correctVelocity(timeStep); comm_.barrier();

    comm_.printf("Max Co = %lf\n", maxCourantNumber(timeStep));

    return 0.;
}

Scalar FractionalStep::maxCourantNumber(Scalar timeStep) const
{
    Scalar maxCo = 0;

    for(const Face &face: grid_.interiorFaces())
    {
        Vector2D sf = face.outwardNorm(face.lCell().centroid());
        Vector2D rc = face.rCell().centroid() - face.lCell().centroid();

        maxCo = std::max(maxCo, fabs(dot(u(face), sf)/dot(rc, sf)));
    }

    return comm_.max(maxCo*timeStep);
}

Scalar FractionalStep::computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
{
    Scalar co = maxCourantNumber(prevTimeStep);
    Scalar lambda1 = 0.1, lambda2 = 1.2;

    return comm_.min(
                std::min(
                    std::min(maxCo/co*prevTimeStep, (1 + lambda1*maxCo/co)*prevTimeStep),
                    std::min(lambda2*prevTimeStep, maxTimeStep_)
                    ));
}

//- Protected methods

Scalar FractionalStep::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rho*u, u) + ib_.eqns(u) == fv::laplacian(mu, u) - fv::source(gradP));
    Scalar error = uEqn_.solve();
    grid_.sendMessages(comm_, u);

    computeFaceVelocities(timeStep);

    return error;
}

Scalar FractionalStep::solvePEqn(Scalar timeStep)
{
    computeMassSource(timeStep);

    pEqn_ = (fv::laplacian(timeStep/rho, dp) + ib_.eqns(dp) == divUStar);
    Scalar error = pEqn_.solveWithGuess();

    grid_.sendMessages(comm_, dp);

    for(const Cell& cell: p.grid.activeCells())
        p(cell) += dp(cell);

    grid_.sendMessages(comm_, p);

    gradP.savePreviousTimeStep(timeStep, 1);

    p.setBoundaryFaces();
    fv::computeInverseWeightedGradient(rho, p, gradP);

    dp.setBoundaryFaces();
    fv::computeInverseWeightedGradient(rho, dp, gradDp);

    return error;
}

void FractionalStep::correctVelocity(Scalar timeStep)
{
    const VectorFiniteVolumeField& gradP0 = gradP.prev(0);

    for(const Cell &cell: grid_.cellZone("fluid"))
        u(cell) -= timeStep/rho(cell)*(gradP(cell) - gradP0(cell));

    for(const Face &face: grid_.interiorFaces())
        u(face) -= timeStep/rho(face)*gradDp(face);

    for(const Face &face: grid_.boundaryFaces())
    {
        switch(u.boundaryType(face))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::SYMMETRY:
        {
            const Vector2D nWall = face.outwardNorm(face.lCell().centroid());
            u(face) = u(face.lCell()) - dot(u(face.lCell()), nWall)*nWall/nWall.magSqr();
        }
            break;

        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            u(face) -= timeStep/rho(face)*gradDp(face);
            break;
        };
    }
    //- Because faces are corrected, no need to communicate the velocity corrections
}

void FractionalStep::computeFaceVelocities(Scalar timeStep)
{
    for(const Face& face: u.grid.interiorFaces())
    {
        const Cell &lCell = face.lCell();
        const Cell &rCell = face.rCell();
        const Scalar g = rCell.volume()/(rCell.volume() + lCell.volume());

        u(face) = g*(u(lCell) + timeStep/rho(lCell)*gradP(lCell))
                + (1. - g)*(u(rCell) + timeStep/rho(rCell)*gradP(rCell))
                - timeStep/rho(face)*gradP(face);
    }

    for(const Face& face: u.grid.boundaryFaces())
    {
        const Cell &cell = face.lCell();

        switch(u.boundaryType(face))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            u(face) = u(cell) + timeStep/rho(cell)*gradP(cell)
                    - timeStep/rho(face)*gradP(face);
            break;

        case VectorFiniteVolumeField::SYMMETRY:
        {
            const Vector2D nWall = face.outwardNorm(cell.centroid());
            u(face) = u(cell) - dot(u(cell), nWall)*nWall/nWall.magSqr();
        }
            break;

            //default:
            //throw Exception("FractionalStep", "computeFaceVelocities", "unrecognized boundary condition type.");
        };
    }
}

void FractionalStep::computeMassSource(Scalar timeStep)
{
    divUStar.fill(0.);

    for(const Cell &cell: grid_.cellZone("fluid"))
    {
        for(const InteriorLink &nb: cell.neighbours())
            divUStar(cell) += dot(u(nb.face()), nb.outwardNorm());

        for(const BoundaryLink &bd: cell.boundaries())
            divUStar(cell) += dot(u(bd.face()), bd.outwardNorm());
    }
}


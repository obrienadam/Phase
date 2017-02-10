#include "FractionalStepExperimental.h"
#include "FiniteVolumeEquation.h"
#include "CrankNicolson.h"
#include "SourceEvaluation.h"
#include "GradientEvaluation.h"
#include "FaceInterpolation.h"

FractionalStepExperimental::FractionalStepExperimental(const Input &input, const Communicator &comm, FiniteVolumeGrid2D &grid)
    :
      Solver(input, comm, grid),
      u(addVectorField(input, "u")),
      gradP(addVectorField(input, "gradP")),
      gradPhi(addVectorField("gradPhi")),
      p(addScalarField(input, "p")),
      phi(addScalarField("phi")),
      rho(addScalarField("rho")),
      mu(addScalarField("mu")),
      divUStar(addScalarField("divUStar")),
      uEqn_(input, comm, u, "uEqn"),
      phiEqn_(input, comm, phi, "phiEqn")
{
    rho.fill(input.caseInput().get<Scalar>("Properties.rho", 1.));
    mu.fill(input.caseInput().get<Scalar>("Properties.mu", 1.));

    rho.savePreviousTimeStep(0., 1);
    mu.savePreviousTimeStep(0., 1);

    phi.copyBoundaryTypes(p);

    //- All active cells to fluid cells
    grid_.createCellZone("fluid", grid_.getCellIds(grid_.localActiveCells()));

    //- Create ib zones if any
    ib_.initCellZones(comm);

    //- Compute the global cell ordering (for lin alg)
    grid_.computeGlobalOrdering(comm_);
}

std::string FractionalStepExperimental::info() const
{
    return Solver::info()
            + "Type: 2nd order fractional-step (experimental)\n"
            + "Advection time-marching: Crank-Nicolson\n"
            + "Diffusion time-marching: Adams-Bashforth\n";
}

Scalar FractionalStepExperimental::solve(Scalar timeStep)
{
    solveUEqn(timeStep);
    solvePhiEqn(timeStep);
    correctVelocity(timeStep);

    comm_.printf("Max Co = %lf\n", maxCourantNumber(timeStep));

    return 0.;
}

Scalar FractionalStepExperimental::maxCourantNumber(Scalar timeStep) const
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

Scalar FractionalStepExperimental::computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
{
    Scalar co = maxCourantNumber(prevTimeStep);
    Scalar lambda1 = 0.1, lambda2 = 1.2;

    return comm_.min(
                std::min(
                    std::min(maxCo/co*prevTimeStep, (1 + lambda1*maxCo/co)*prevTimeStep),
                    std::min(lambda2*prevTimeStep, maxTimeStep_)
                    ));
}

Scalar FractionalStepExperimental::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    uEqn_ = (fv::ddt(rho, u, timeStep) + cn::div(rho, u, u, 0.) + ib_.eqns(u) == cn::laplacian(mu, u, 0.) - fv::source(gradP));

    Scalar error = uEqn_.solve();
    grid_.sendMessages(comm_, u);

    computeFaceVelocities(timeStep);

    return error;
}

Scalar FractionalStepExperimental::solvePhiEqn(Scalar timeStep)
{
    phi.savePreviousTimeStep(timeStep, 1);

    computeMassSource(timeStep);

    phiEqn_ = (fv::laplacian(1., phi) + ib_.eqns(phi) == divUStar);
    Scalar error = phiEqn_.solve();

    grid_.sendMessages(comm_, phi);

    phi.setBoundaryFaces();
    fv::computeInverseWeightedGradient(rho, phi, gradPhi);

    for(const Cell& cell: grid_.localActiveCells())
        gradP(cell) += gradPhi(cell)*rho(cell)/timeStep;

    return error;
}

void FractionalStepExperimental::correctVelocity(Scalar timeStep)
{
    for(const Cell& cell: grid_.localActiveCells())
        u(cell) -= gradPhi(cell);

    for(const Face& face: grid_.interiorFaces())
        u(face) -= gradPhi(face);
}

void FractionalStepExperimental::computeFaceVelocities(Scalar timeStep)
{
    u.setBoundaryFaces();
    fv::interpolateFaces(fv::INVERSE_DISTANCE, u);
    return;

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

void FractionalStepExperimental::computeMassSource(Scalar timeStep)
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

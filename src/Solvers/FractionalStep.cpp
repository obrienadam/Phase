#include "FractionalStep.h"
#include "CrankNicolson.h"
#include "AdamsBashforth.h"
#include "GradientEvaluation.h"
#include "FaceInterpolation.h"

FractionalStep::FractionalStep(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Solver(grid, input),
      u(addVectorField(input, "u")),
      sg(addVectorField("sg")),
      gradP(addVectorField("gradP")),
      p(addScalarField(input, "p")),
      rho(addScalarField("rho")),
      mu(addScalarField("mu")),
      divUStar(addScalarField("uStar")),
      uEqn_(u, "uEqn", SparseMatrix::IncompleteLUT),
      pEqn_(p, "pEqn", SparseMatrix::IncompleteLUT)
{
    rho.fill(input.caseInput().get<Scalar>("Properties.rho", 1.));
    mu.fill(input.caseInput().get<Scalar>("Properties.mu", 1.));
    g_ = Vector2D(input.caseInput().get<std::string>("Properties.g"));

    volumeIntegrators_ = VolumeIntegrator::initVolumeIntegrators(input, *this);
    forceIntegrators_ = ForceIntegrator::initForceIntegrators(input, p, rho, mu, u);

    p.initNodes();

    uEqn_.matrix().setFill(1);
    pEqn_.matrix().setFill(3);

    interpolateFaces(fv::INVERSE_VOLUME, u);
}

std::string FractionalStep::info() const
{
    return Solver::info()
            + "Type: 2nd order fractional-step\n"
            + "Advection time-marching: Crank-Nicolson\n"
            + "Diffusion time-marching: Adams-Bashforth\n"
            + "Gravity: " + g_.toString() + "\n";
}

Scalar FractionalStep::solve(Scalar timeStep)
{
    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);

    for(const ForceIntegrator &fi: forceIntegrators_)
        fi.integrate();

    for(const VolumeIntegrator &vi: volumeIntegrators_)
        vi.integrate();

    printf("Max Co = %lf\n", courantNumber(timeStep));

    return 0.;
}

Scalar FractionalStep::computeMaxTimeStep(Scalar maxCo) const
{
    Scalar maxTimeStepSqr = std::numeric_limits<Scalar>::infinity(), maxCoSqr = maxCo*maxCo;

    for(const Cell &cell: u.grid.fluidCells())
        for(const InteriorLink &nb: cell.neighbours())
        {
            Scalar deltaXSqr = nb.rCellVec().magSqr();
            Scalar magUSqr = u(nb.face()).magSqr();

            maxTimeStepSqr = std::min(maxTimeStepSqr, maxCoSqr*deltaXSqr/magUSqr);
        }

    return sqrt(maxTimeStepSqr);
}

//- Protected methods

Scalar FractionalStep::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    sg = fv::gravity(rho, g_);

    uEqn_ = (fv::ddt(rho, u, timeStep) + cn::div(rho*u, u) + ib_.eqns(u) == ab::laplacian(mu, u) + fv::source(sg - gradP));

    Scalar error = uEqn_.solve();

    for(const Cell& cell: u.grid.fluidCells())
        u(cell) += timeStep*(gradP(cell) - sg(cell))/rho(cell);

    interpolateFaces(fv::INVERSE_VOLUME, u);

    return error;
}

Scalar FractionalStep::solvePEqn(Scalar timeStep)
{
    computeMassSource(timeStep);

    pEqn_ = (fv::laplacian(p) + ib_.eqns(p) == divUStar);
    Scalar error = pEqn_.solve(p.sparseVector());

    computeGradient(fv::GREEN_GAUSS_CELL_CENTERED, p, gradP);

    return error;
}

void FractionalStep::correctVelocity(Scalar timeStep)
{
    for(const Cell &cell: grid_.fluidCells())
        u(cell) += timeStep/rho(cell)*(sg(cell) - gradP(cell));

    for(const Face &face: grid_.interiorFaces())
    {
        const Cell &lCell = face.lCell();
        const Cell &rCell = face.rCell();
        Vector2D rc = rCell.centroid() - lCell.centroid();

        u(face) += timeStep/rho(face)*(p(lCell) - p(rCell))*rc/dot(rc, rc) + timeStep*g_;
    }

    for(const Face &face: grid_.boundaryFaces())
    {
        const Cell &cell = face.lCell();
        Vector2D rf = face.centroid() - cell.centroid();

        if(u.boundaryType(face.id()) == VectorFiniteVolumeField::FIXED)
            continue;
        else if(u.boundaryType(face.id()) == VectorFiniteVolumeField::SYMMETRY) // This is not correct
            continue;
        else
            u(face) += timeStep/rho(face)*(p(cell) - p(face))*rf/dot(rf, rf) + timeStep*g_;
    }
}

void FractionalStep::computeMassSource(Scalar timeStep)
{
    divUStar.fill(0.);

    for(const Cell &cell: grid_.fluidCells())
    {
        for(const InteriorLink &nb: cell.neighbours())
        {
            const Vector2D& sf = nb.outwardNorm();
            divUStar(cell) += rho(nb.face())/timeStep*dot(u(nb.face()), sf);
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Vector2D& sf = bd.outwardNorm();
            divUStar(cell) += rho(bd.face())/timeStep*dot(u(bd.face()), sf);
        }
    }
}

Scalar FractionalStep::courantNumber(Scalar timeStep)
{
    Scalar maxCoSqr = 0., timeStepSqr = timeStep*timeStep;

    for(const Cell &cell: u.grid.fluidCells())
        for(const InteriorLink &nb: cell.neighbours())
        {
            Scalar deltaXSqr = nb.rCellVec().magSqr();
            Scalar magUSqr = u(nb.face()).magSqr();

            maxCoSqr = std::max(maxCoSqr, magUSqr*timeStepSqr/deltaXSqr);
        }

    return sqrt(maxCoSqr);
}


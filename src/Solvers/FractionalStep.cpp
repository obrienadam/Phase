#include "FractionalStep.h"
#include "CrankNicolson.h"
#include "AdamsBashforth.h"

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

    uEqn_.matrix().setFill(2);
    pEqn_.matrix().setFill(3);

    u.savePreviousTimeStep(0., 1);
    computeAdvectingVelocity(0.);
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

    printf("Current time step = %lf\n", timeStep);
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
            Scalar magUSqr = u.faces()[nb.face().id()].magSqr();

            maxTimeStepSqr = std::min(maxTimeStepSqr, maxCoSqr*deltaXSqr/magUSqr);
        }

    return sqrt(maxTimeStepSqr);
}

//- Protected methods

Scalar FractionalStep::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    uEqn_ = (fv::ddt(rho, u, timeStep) + cn::div(rho*u, u) + ib_.eqns(u) == ab::laplacian(mu, u) - fv::source(gradP));

    Scalar error = uEqn_.solve();
    computeAdvectingVelocity(timeStep);

    return error;
}

Scalar FractionalStep::solvePEqn(Scalar timeStep)
{
    p.savePreviousTimeStep(timeStep, 1);
    divUStar.fill(0.);

    for(const Cell &cell: grid_.fluidCells())
    {
        Scalar rhof, p1;

        const Scalar p0 = p[cell.id()];

        for(const InteriorLink &nb: cell.neighbours())
        {
            rhof = rho.faces()[nb.face().id()];

            const Vector2D& uf = u.faces()[nb.face().id()];
            const Vector2D& sf = nb.outwardNorm();
            const Vector2D& rc = nb.rCellVec();

            p1 = p[nb.cell().id()];

            divUStar[cell.id()] += rhof/timeStep*dot(uf, sf) + (p1 - p0)*dot(rc, sf)/dot(rc, rc);
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            rhof = rho.faces()[bd.face().id()];

            const Vector2D& uf = u.faces()[bd.face().id()];
            const Vector2D& sf = bd.outwardNorm();
            const Vector2D& rf = bd.rFaceVec();

            p1 = p.faces()[bd.face().id()];

            divUStar[cell.id()] += rhof/timeStep*dot(uf, sf) + (p1 - p0)*dot(rf, sf)/dot(rf, rf);
        }
    }

    pEqn_ = (fv::laplacian(p) + ib_.eqns(p) == divUStar);
    Scalar error = pEqn_.solve();

    gradP.savePreviousTimeStep(timeStep, 1);

    interpolateFaces(p);
    gradP = grad(p);

    return error;
}

void FractionalStep::correctVelocity(Scalar timeStep)
{
    const ScalarFiniteVolumeField &p0 = p.prev(0);
    const VectorFiniteVolumeField &gradP0 = gradP.prev(0);


    for(const Cell &cell: grid_.fluidCells())
    {
        u[cell.id()] -= timeStep/rho[cell.id()]*gradP[cell.id()] - timeStep/rho[cell.id()]*gradP0[cell.id()];
    }

    for(const Face &face: grid_.interiorFaces())
    {
        const Cell &lCell = face.lCell();
        const Cell &rCell = face.rCell();
        Vector2D rc = rCell.centroid() - lCell.centroid();

        u.faces()[face.id()] -= timeStep/rho.faces()[face.id()]*(p[rCell.id()] - p[lCell.id()]
                - p0[rCell.id()] + p0[lCell.id()])*rc/rc.magSqr();
    }

    for(const Face &face: grid_.boundaryFaces())
    {
        const Cell &cell = face.lCell();
        Vector2D rf = face.centroid() - cell.centroid();

        if(u.boundaryType(face.id()) == VectorFiniteVolumeField::FIXED)
            continue;
        else if(u.boundaryType(face.id()) == VectorFiniteVolumeField::SYMMETRY)
            continue;
        else
            u.faces()[face.id()] -= timeStep/rho.faces()[face.id()]*(p.faces()[face.id()] - p[cell.id()]
                    - p0.faces()[face.id()] + p0[cell.id()])*rf/rf.magSqr();
    }
}

void FractionalStep::computeAdvectingVelocity(Scalar timeStep)
{
    for(const Face &face: u.grid.interiorFaces())
    {
        const Cell &lCell = face.lCell();
        const Cell &rCell = face.rCell();

        Scalar g = rCell.volume()/(lCell.volume() + rCell.volume());
        Vector2D rc = rCell.centroid() - lCell.centroid();

        u.faces()[face.id()] = g*(u[lCell.id()] + timeStep*gradP[lCell.id()]/rho[lCell.id()])
                + (1. - g)*(u[rCell.id()] + timeStep*gradP[rCell.id()]/rho[rCell.id()])
                - timeStep/rho.faces()[face.id()]*(p[rCell.id()] - p[lCell.id()])*rc/rc.magSqr();
    }

    for(const Face &face: u.grid.boundaryFaces())
    {
        const Cell &cell = face.lCell();

        switch(u.boundaryType(face.id()))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::NORMAL_GRADIENT: case VectorFiniteVolumeField::OUTFLOW:
        {
            Vector2D rf = face.centroid() - cell.centroid();

            u.faces()[face.id()] = u[cell.id()] + timeStep*gradP[cell.id()]/rho[cell.id()]
                    - timeStep/rho.faces()[face.id()]*(p.faces()[face.id()] - p[cell.id()])*rf/rf.magSqr();
        }
            break;

        case VectorFiniteVolumeField::SYMMETRY:
        {
            Vector2D nWall = face.outwardNorm(face.lCell().centroid());
            u.faces()[face.id()] = u[face.lCell().id()] - dot(u[face.lCell().id()], nWall)*nWall/nWall.magSqr();
            break;
        }

        default:
            throw Exception("FractionalStep", "computeAdvectingVelocity", "unrecongnized boundary condition type.");
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
            Scalar magUSqr = u.faces()[nb.face().id()].magSqr();

            maxCoSqr = std::max(maxCoSqr, magUSqr*timeStepSqr/deltaXSqr);
        }

    return sqrt(maxCoSqr);
}


#include "FractionalStep.h"
#include "CrankNicolson.h"
#include "AdamsBashforth.h"

FractionalStep::FractionalStep(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Solver(grid, input),
      u(addVectorField(input, "u")),
      p(addScalarField(input, "p")),
      rho(addScalarField("rho")),
      mu(addScalarField("mu")),
      divUStar(addScalarField("uStar")),
      uEqn_(u, "uEqn", SparseMatrix::IncompleteLUT),
      pEqn_(p, "pEqn", SparseMatrix::IncompleteLUT),
      ib_(input, *this)
{
    rho.fill(input.caseInput().get<Scalar>("Properties.rho", 1.));
    mu.fill(input.caseInput().get<Scalar>("Properties.mu", 1.));

    uEqn_.matrix().setFill(2);
    pEqn_.matrix().setFill(3);

    u.savePreviousTimeStep(0., 1);
}

Scalar FractionalStep::solve(Scalar timeStep)
{
    u.savePreviousTimeStep(0., 1);

    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);

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
    uEqn_ = (fv::ddt(rho, u, timeStep) + cn::div(rho*u, u) + ib_.eqns(u) == ab::laplacian(mu, u));
    Scalar error = uEqn_.solve();
    interpolateFaces(u);
    return error;
}

Scalar FractionalStep::solvePEqn(Scalar timeStep)
{
    divUStar.fill(0.);

    for(const Cell &cell: grid_.fluidCells())
    {
        for(const InteriorLink &nb: cell.neighbours())
            divUStar[cell.id()] += rho[cell.id()]/timeStep*dot(u.faces()[nb.face().id()], nb.outwardNorm());

        for(const BoundaryLink &bd: cell.boundaries())
            divUStar[cell.id()] += rho[cell.id()]/timeStep*dot(u.faces()[bd.face().id()], bd.outwardNorm());
    }

    pEqn_ = (fv::laplacian(p) + ib_.eqns(p) == divUStar);
    Scalar error = pEqn_.solve();

    interpolateFaces(p);

    return error;
}

void FractionalStep::correctVelocity(Scalar timeStep)
{
    VectorFiniteVolumeField gradP = grad(p);

    for(const Cell &cell: grid_.fluidCells())
        u[cell.id()] -= timeStep/rho[cell.id()]*gradP[cell.id()];

    for(const Face &face: grid_.interiorFaces())
    {
        const Cell &lCell = face.lCell();
        const Cell &rCell = face.rCell();
        Vector2D rc = rCell.centroid() - lCell.centroid();

        u.faces()[face.id()] -= timeStep/rho.faces()[face.id()]*(p[rCell.id()] - p[lCell.id()])*rc/dot(rc, rc);
    }

    for(const Face &face: grid_.boundaryFaces())
    {
        const Cell &cell = face.lCell();
        Vector2D rf = face.centroid() - cell.centroid();

        u.faces()[face.id()] -= timeStep/rho.faces()[face.id()]*(p.faces()[face.id()] - p[cell.id()])*rf/dot(rf, rf);
    }
}

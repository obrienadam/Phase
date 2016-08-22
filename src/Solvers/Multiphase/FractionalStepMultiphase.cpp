#include "FractionalStepMultiphase.h"
#include "Cicsam.h"
#include "CrankNicolson.h"
#include "AdamsBashforth.h"

FractionalStepMultiphase::FractionalStepMultiphase(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      FractionalStep(grid, input),
      gamma(addScalarField(input, "gamma")),
      gradGamma(addVectorField("gradGamma")),
      ft(addVectorField("ft")),
      gammaEqn_(gamma, "gammaEqn"),
      surfaceTensionForce_(input, *this)
{
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1");
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2");
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1");
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2");

    setInitialConditions(input);

    surfaceTensionForce_.compute();
    computeRho();
    computeMu();
}

Scalar FractionalStepMultiphase::solve(Scalar timeStep)
{
    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);
    solveGammaEqn(timeStep);

    printf("Max Co = %lf\n", courantNumber(timeStep));

    return 0.;
}

//- Protected methods

Scalar FractionalStepMultiphase::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    ft = surfaceTensionForce_.compute();
    sg = fv::gravity(rho, g_);

    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rho*u, u) + ib_.eqns(u)
             == fv::laplacian(mu, u) - fv::source(gradP) + fv::source(sg) + fv::source(ft));

    Scalar error = uEqn_.solve();
    computeAdvectingVelocity(timeStep);

    return error;
}

Scalar FractionalStepMultiphase::solveGammaEqn(Scalar timeStep)
{
    gamma.savePreviousTimeStep(timeStep, 1);

    interpolateFaces(gamma);

    //gammaEqn_ = (cicsam::div(u, gamma, timeStep, cicsam::HC) + ib_.eqns(gamma) == 0.);
    gammaEqn_ = (cicsam::div(u, gradGamma, ib_.ibObjs(), gamma, timeStep, cicsam::HC) == 0.);

    Scalar error = gammaEqn_.solve();

    computeRho();
    computeMu();

    return error;
}

void FractionalStepMultiphase::computeMassSource(Scalar timeStep)
{
    p.savePreviousTimeStep(timeStep, 1);
    divUStar.fill(0.);
    const Scalar sigma = surfaceTensionForce_.sigma();

    for(const Cell &cell: grid_.fluidCells())
    {
        Scalar rhof, kf, p1, g1;

        const Scalar p0 = p[cell.id()];
        const Scalar g0 = gamma[cell.id()];

        for(const InteriorLink &nb: cell.neighbours())
        {
            rhof = rho.faces()[nb.face().id()];
            kf = surfaceTensionForce_.kappa().faces()[nb.face().id()];

            const Vector2D& uf = u.faces()[nb.face().id()];
            const Vector2D& sf = nb.outwardNorm();
            const Vector2D& rc = nb.rCellVec();

            p1 = p[nb.cell().id()];
            g1 = gamma[nb.cell().id()];

            divUStar[cell.id()] += rhof/timeStep*dot(uf, sf) + (p1 - p0)*dot(rc, sf)/dot(rc, rc);// - rhof*dot(g_, sf) - sigma*kf*(g1 - g0)*dot(rc, sf)/dot(rc, rc);
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            rhof = rho.faces()[bd.face().id()];
            kf = surfaceTensionForce_.kappa().faces()[bd.face().id()];

            const Vector2D& uf = u.faces()[bd.face().id()];
            const Vector2D& sf = bd.outwardNorm();
            const Vector2D& rf = bd.rFaceVec();

            p1 = p.faces()[bd.face().id()];
            g1 = gamma.faces()[bd.face().id()];

            divUStar[cell.id()] += rhof/timeStep*dot(uf, sf) + (p1 - p0)*dot(rf, sf)/dot(rf, rf);// - rhof*dot(g_, sf) - sigma*kf*(g1 - g0)*dot(rf, sf)/dot(rf, rf);
        }
    }
}

void FractionalStepMultiphase::computeAdvectingVelocity(Scalar timeStep)
{
    //FractionalStep::computeAdvectingVelocity(timeStep);
    const Scalar sigma = surfaceTensionForce_.sigma();

    for(const Face &face: u.grid.interiorFaces())
    {
        const Cell &lCell = face.lCell();
        const Cell &rCell = face.rCell();

        Scalar g = rCell.volume()/(lCell.volume() + rCell.volume());
        Vector2D rc = rCell.centroid() - lCell.centroid();

        const Scalar rhof = rho.faces()[face.id()];
        const Scalar kf = surfaceTensionForce_.kappa().faces()[face.id()];
        const Scalar p0 = p[lCell.id()], p1 = p[rCell.id()];
        const Scalar g0 = gamma[lCell.id()], g1 = gamma[rCell.id()];

        u.faces()[face.id()] = g*(u[lCell.id()] + timeStep*(gradP[lCell.id()] - sg[lCell.id()] - ft[lCell.id()])/rho[lCell.id()])
                + (1. - g)*(u[rCell.id()] + timeStep*(gradP[rCell.id()] - sg[rCell.id()] - ft[rCell.id()])/rho[rCell.id()])
                - timeStep/rhof*(p1 - p0)*rc/rc.magSqr()
                + timeStep/rhof*g_;
                + timeStep/rhof*sigma*kf*(g1 - g0)*rc/rc.magSqr();
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
            const Scalar rhof = rho.faces()[face.id()];
            const Scalar kf = surfaceTensionForce_.kappa().faces()[face.id()];
            const Scalar p0 = p[cell.id()], p1 = p.faces()[face.id()];
            const Scalar g0 = gamma[cell.id()], g1 = gamma.faces()[face.id()];

            u.faces()[face.id()] = u[cell.id()]
                    + timeStep*(gradP[cell.id()] - sg[cell.id()] - ft[cell.id()])/rho[cell.id()]
                    - timeStep/rhof*(p1 - p0)*rf/rf.magSqr()
                    + timeStep/rhof*g_
                    + timeStep/rhof*sigma*kf*(g1 - g0)*rf/rf.magSqr();
        }
            break;

        case VectorFiniteVolumeField::SYMMETRY:
        {
            Vector2D nWall = face.outwardNorm(face.lCell().centroid());
            u.faces()[face.id()] = u[face.lCell().id()] - dot(u[face.lCell().id()], nWall)*nWall/nWall.magSqr();
            break;
        }

        default:
            throw Exception("FractionalStepMultiphase", "computeAdvectingVelocity", "unrecongnized boundary condition type.");
        }
    }
}

void FractionalStepMultiphase::computeRho()
{
    const ScalarFiniteVolumeField &w = surfaceTensionForce_.gammaTilde();

    for(const Cell &cell: grid_.activeCells())
        rho[cell.id()] = (1. - w[cell.id()])*rho1_ + w[cell.id()]*rho2_;

    interpolateFaces(rho);
    //harmonicInterpolateFaces(rho);
}

void FractionalStepMultiphase::computeMu()
{
    const ScalarFiniteVolumeField &w = surfaceTensionForce_.gammaTilde();

    for(const Cell &cell: grid_.activeCells())
        mu[cell.id()] = (1. - w[cell.id()])*mu1_ + w[cell.id()]*mu2_;

    interpolateFaces(mu);
    //harmonicInterpolateFaces(mu);
}

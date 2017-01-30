#include "Piso.h"
#include "Exception.h"
#include "CrankNicolson.h"
#include "FaceInterpolation.h"
#include "GradientEvaluation.h"
#include "SourceEvaluation.h"

Piso::Piso(const Input &input, const Communicator &comm, FiniteVolumeGrid2D& grid)
    :
      Solver(input, comm, grid),
      u(addVectorField(input, "u")),
      gradP(addVectorField("gradP")),
      gradPCorr(addVectorField("gradPCorr")),
      p(addScalarField(input, "p")),
      pCorr(addScalarField("pCorr")),
      rho(addScalarField("rho")),
      mu(addScalarField("mu")),
      m(addScalarField("m")),
      d(addScalarField("d")),
      uEqn_(input, comm, u, "uEqn"),
      pCorrEqn_(input, comm, pCorr, "pCorrEqn")
{
    rho.fill(input.caseInput().get<Scalar>("Properties.rho", 1.));
    mu.fill(input.caseInput().get<Scalar>("Properties.mu", 1.));

    rho.savePreviousTimeStep(0, 1);
    mu.savePreviousTimeStep(0, 1);

    nInnerIterations_ = input.caseInput().get<size_t>("Solver.numInnerIterations");
    nPCorrections_ = input.caseInput().get<size_t>("Solver.numPressureCorrections");
    momentumOmega_ = input.caseInput().get<Scalar>("Solver.momentumRelaxation");
    pCorrOmega_ = input.caseInput().get<Scalar>("Solver.pressureCorrectionRelaxation");
    pCorr.copyBoundaryTypes(p);

    //- All active cells to fluid cells
    grid_.createCellZone("fluid", grid_.getCellIds(grid_.activeCells()));
    ib_.initCellZones();
    grid_.computeOrdering();

    volumeIntegrators_ = VolumeIntegrator::initVolumeIntegrators(input, *this);
    forceIntegrators_ = ForceIntegrator::initForceIntegrators(input, p, rho, mu, u);
}

Scalar Piso::solve(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    for(size_t innerIter = 0; innerIter < nInnerIterations_; ++innerIter)
    {
        u.savePreviousIteration();
        solveUEqn(timeStep);

        for(size_t pCorrIter = 0; pCorrIter < nPCorrections_; ++pCorrIter)
        {
            solvePCorrEqn();
            correctPressure();
            correctVelocity();
        }
    }

    printf("Max Co = %lf\n", maxCourantNumber(timeStep));

    return 0.;
}

Scalar Piso::maxCourantNumber(Scalar timeStep) const
{
    Scalar maxCo = 0;

    for(const Face &face: grid_.interiorFaces())
    {
        Vector2D sf = face.outwardNorm(face.lCell().centroid());
        Vector2D rc = face.rCell().centroid() - face.lCell().centroid();

        maxCo = std::max(maxCo, fabs(dot(u(face), sf)/dot(rc, sf)));
    }

    return maxCo*timeStep;
}

Scalar Piso::computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
{
    Scalar co = maxCourantNumber(prevTimeStep);
    Scalar lambda1 = 0.05, lambda2 = 1.1;

    return std::min(
                std::min(maxCo/co*prevTimeStep, (1 + lambda1*maxCo/co)*prevTimeStep),
                std::min(lambda2*prevTimeStep, maxTimeStep_)
                );
}

//- Protected methods

Scalar Piso::solveUEqn(Scalar timeStep)
{
    uEqn_ = (fv::ddt(rho, u, timeStep) + cn::div(rho, u, u, 1.5) + ib_.eqns(u) == cn::laplacian(mu, u, 0.5) - fv::source(gradP));
    uEqn_.relax(momentumOmega_);

    Scalar error = uEqn_.solve();

    rhieChowInterpolation();

    return error;
}

Scalar Piso::solvePCorrEqn()
{
    for(const Cell& cell: m.grid.cellZone("fluid"))
    {
        m(cell) = 0.;

        for(const InteriorLink& nb: cell.neighbours())
            m(cell) += dot(u(nb.face()), nb.outwardNorm());

        for(const BoundaryLink& bd: cell.boundaries())
            m(cell) += dot(u(bd.face()), bd.outwardNorm());
    }

    pCorrEqn_ = (fv::laplacian(d, pCorr) + ib_.eqns(pCorr) == m);
    Scalar error = pCorrEqn_.solve();

    pCorr.setBoundaryFaces();
    fv::computeInverseWeightedGradient(rho, pCorr, gradPCorr);

    return error;
}

void Piso::rhieChowInterpolation()
{   
    VectorFiniteVolumeField& uStar = u.prevIter();
    VectorFiniteVolumeField& uPrev = u.prev(0);
    ScalarFiniteVolumeField& rhoPrev = rho.prev(0);
    const Scalar dt = u.prevTimeStep(0);

    d.fill(0.);
    for(const Cell& cell: d.grid.cellZone("fluid"))
    {
        Vector2D coeff = uEqn_.get(cell, cell);
        d(cell) = cell.volume()/(0.5*(coeff.x + coeff.y));
    }

    interpolateFaces(fv::INVERSE_VOLUME, d);

    for(const Face& face: u.grid.interiorFaces())
    {
        const Cell& cellP = face.lCell();
        const Cell& cellQ = face.rCell();

        const Scalar dP = d(cellP);
        const Scalar dQ = d(cellQ);
        const Scalar df = d(face);

        const Scalar rhoP0 = rhoPrev(cellP);
        const Scalar rhoQ0 = rhoPrev(cellQ);
        const Scalar rhof0 = rhoPrev(face);

        const Scalar g = cellQ.volume()/(cellP.volume() + cellQ.volume());


        u(face) = g*u(cellP) + (1. - g)*u(cellQ)
                + (1. - momentumOmega_)*(uStar(face) - (g*uStar(cellP) + (1. - g)*uStar(cellQ)))
                + (rhof0*df*uPrev(face) - (g*rhoP0*dP*uPrev(cellP) + (1. - g)*rhoQ0*dQ*uPrev(cellQ)))/dt //- This term is very important!
                - df*gradP(face) + (g*dP*gradP(cellP) + (1. - g)*dQ*gradP(cellQ));
    }

    for(const Face& face: u.grid.boundaryFaces())
    {
        const Scalar df = d(face);
        const Scalar rhoP0 = rhoPrev(face.lCell());
        const Scalar rhof0 = rhoPrev(face);

        switch(u.boundaryType(face))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            u(face) = u(face.lCell())
                    + rhof0*df*uPrev(face) - rhoP0*d(face.lCell())*uPrev(face.lCell())
                    - df*gradP(face) + d(face.lCell())*gradP(face.lCell());
            break;

        case VectorFiniteVolumeField::SYMMETRY:
        {
            const Vector2D nWall = face.outwardNorm(face.lCell().centroid());

            u(face) = u(face.lCell()) - dot(u(face.lCell()), nWall)*nWall/nWall.magSqr();
        }
            break;

        default:
            throw Exception("Piso", "rhieChowInterpolation", "unrecognized boundary condition type.");
        }
    }
}

void Piso::correctPressure()
{
    for(const Cell& cell: p.grid.activeCells())
        p(cell) += pCorrOmega_*pCorr(cell);

    p.setBoundaryFaces();
    fv::computeInverseWeightedGradient(rho, p, gradP);
}

void Piso::correctVelocity()
{
    for(const Cell& cell: u.grid.cellZone("fluid"))
        u(cell) -= d(cell)*gradPCorr(cell);

    for(const Face& face: u.grid.interiorFaces())
        u(face) -= d(face)*gradPCorr(face);

    for(const Face& face: u.grid.boundaryFaces())
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
            u(face) -= d(face)*gradPCorr(face);
            break;
        }
    }

    //- Update mass source for the purpose of error checking

    for(const Cell& cell: m.grid.cellZone("fluid"))
    {
        m(cell) = 0.;

        for(const InteriorLink& nb: cell.neighbours())
            m(cell) += dot(u(nb.face()), nb.outwardNorm());

        for(const BoundaryLink& bd: cell.boundaries())
            m(cell) += dot(u(bd.face()), bd.outwardNorm());
    }
}

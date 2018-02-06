#include "Piso.h"
#include "ScalarGradient.h"
#include "Source.h"

Piso::Piso(const Input &input,
           std::shared_ptr<FiniteVolumeGrid2D>& grid)
        :
        Solver(input, grid),
        u(addVectorField(input, "u")),
        p(addScalarField(input, "p")),
        pCorr(addScalarField("pCorr")),
        gradP(addVectorField(std::make_shared<ScalarGradient>(p))),
        gradPCorr(addVectorField(std::make_shared<ScalarGradient>(pCorr))),
        rho(addScalarField("rho")),
        mu(addScalarField("mu")),
        m(addScalarField("m")),
        d(addScalarField("d")),
        uEqn_(input, u, "uEqn"),
        pCorrEqn_(input, pCorr, "pCorrEqn"),
        fluid_(grid->createCellZone("fluid"))
{
    rho.fill(input.caseInput().get<Scalar>("Properties.rho", 1.));
    mu.fill(input.caseInput().get<Scalar>("Properties.mu", 1.));

    rho.savePreviousTimeStep(0, 1);
    mu.savePreviousTimeStep(0, 1);

    //- Solver parameters
    nInnerIterations_ = input.caseInput().get<size_t>("Solver.numInnerIterations");
    nPCorrections_ = input.caseInput().get<size_t>("Solver.numPressureCorrections");
    momentumOmega_ = input.caseInput().get<Scalar>("Solver.momentumRelaxation");
    pCorrOmega_ = input.caseInput().get<Scalar>("Solver.pressureCorrectionRelaxation");

    //- Copy relevant boundary types
    pCorr.copyBoundaryTypes(p);
    ib_->copyBoundaryTypes(p, pCorr);

    //- All active cells to fluid cells
    fluid_.add(grid_->localActiveCells());

    //- Create ib zones if any
    ib_->initCellZones(fluid_);
}

Scalar Piso::solve(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    for (size_t innerIter = 0; innerIter < nInnerIterations_; ++innerIter)
    {
        u.savePreviousIteration();
        solveUEqn(timeStep);

        for (size_t pCorrIter = 0; pCorrIter < nPCorrections_; ++pCorrIter)
        {
            solvePCorrEqn();
            correctVelocity();
        }
    }

    printf("Max Co = %lf\n", maxCourantNumber(timeStep));

    return 0.;
}

Scalar Piso::maxCourantNumber(Scalar timeStep) const
{
    Scalar maxCo = 0;

    for(const Cell& cell: fluid_)
    {
        Scalar co = 0.;

        for(const InteriorLink& nb: cell.neighbours())
            co += std::max(dot(u(nb.face()), nb.outwardNorm()), 0.);

        for(const BoundaryLink& bd: cell.boundaries())
            co += std::max(dot(u(bd.face()), bd.outwardNorm()), 0.);

        co *= timeStep / cell.volume();
        maxCo = std::max(co, maxCo);
    }

    return grid_->comm().max(maxCo);
}

Scalar Piso::computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
{
    Scalar co = maxCourantNumber(prevTimeStep);
    Scalar lambda1 = 0.05, lambda2 = 1.1;

    return grid_->comm().min(
            std::min(
                    std::min(maxCo / co * prevTimeStep, (1 + lambda1 * maxCo / co) * prevTimeStep),
                    std::min(lambda2 * prevTimeStep, maxTimeStep_)
            ));
}

//- Protected methods

Scalar Piso::solveUEqn(Scalar timeStep)
{
    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rho*u, u) + ib_->bcs(u) ==
             fv::laplacian(mu, u) - src::src(gradP, fluid_));

    uEqn_.relax(momentumOmega_);

    Scalar error = uEqn_.solve();

    grid_->sendMessages(u);

    rhieChowInterpolation();

    return error;
}

Scalar Piso::solvePCorrEqn()
{
    m.fill(0);

    for (const Cell &cell: fluid_)
    {
        for (const InteriorLink &nb: cell.neighbours())
            m(cell) += dot(u(nb.face()), nb.outwardNorm());

        for (const BoundaryLink &bd: cell.boundaries())
            m(cell) += dot(u(bd.face()), bd.outwardNorm());
    }

    pCorrEqn_ = (fv::laplacian(d, pCorr) + ib_->bcs(pCorr) == m);

    Scalar error = pCorrEqn_.solve();
    grid_->sendMessages(pCorr);

    pCorr.setBoundaryFaces();
    gradPCorr.compute(fluid_);

    for(const Cell &cell: grid_->localActiveCells())
        p(cell) += pCorrOmega_*pCorr(cell);

    grid_->sendMessages(p);

    p.setBoundaryFaces();
    gradP.compute(fluid_);

    grid_->sendMessages(gradP);

    return error;
}

void Piso::rhieChowInterpolation()
{
    const VectorFiniteVolumeField &uStar = u.prevIteration();
    const VectorFiniteVolumeField &uPrev = u.oldField(0);
    const ScalarFiniteVolumeField &rhoPrev = rho.oldField(0);
    Scalar dt = u.oldTimeStep(0);

    d.fill(0.);
    for (const Cell &cell: grid_->cellZone("fluid"))
    {
        Vector2D coeff = uEqn_.get(cell, cell);
        d(cell) = cell.volume() / (0.5 * (coeff.x + coeff.y));
    }

    grid_->sendMessages(d);
    d.interpolateFaces();

    for (const Face &face: grid_->interiorFaces())
    {
        const Cell &cellP = face.lCell();
        const Cell &cellQ = face.rCell();

        Scalar dP = d(cellP);
        Scalar dQ = d(cellQ);
        Scalar df = d(face);

        Scalar rhoP0 = rhoPrev(cellP);
        Scalar rhoQ0 = rhoPrev(cellQ);
        Scalar rhof0 = rhoPrev(face);

        Scalar g = cellQ.volume() / (cellP.volume() + cellQ.volume());


        u(face) = g * u(cellP) + (1. - g) * u(cellQ)
                  + (1. - momentumOmega_) * (uStar(face) - (g * uStar(cellP) + (1. - g) * uStar(cellQ)))
                  +
                  (rhof0 * df * uPrev(face) - (g * rhoP0 * dP * uPrev(cellP) + (1. - g) * rhoQ0 * dQ * uPrev(cellQ))) /
                  dt //- This term is very important!
                  - df * gradP(face) + (g * dP * gradP(cellP) + (1. - g) * dQ * gradP(cellQ));
    }

    for (const Face &face: u.grid()->boundaryFaces())
    {
        Scalar df = d(face);
        Scalar rhoP0 = rhoPrev(face.lCell());
        Scalar rhof0 = rhoPrev(face);

        switch (u.boundaryType(face))
        {
            case VectorFiniteVolumeField::FIXED:
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                u(face) = u(face.lCell())
                          + rhof0 * df * uPrev(face) - rhoP0 * d(face.lCell()) * uPrev(face.lCell())
                          - df * gradP(face) + d(face.lCell()) * gradP(face.lCell());
                break;

            case VectorFiniteVolumeField::SYMMETRY:
            {
                const Vector2D nWall = face.outwardNorm(face.lCell().centroid());

                u(face) = u(face.lCell()) - dot(u(face.lCell()), nWall) * nWall / nWall.magSqr();
            }
                break;

                //default:
                //    throw Exception("Piso", "rhieChowInterpolation", "unrecognized boundary condition type.");
        }
    }
}

void Piso::correctVelocity()
{
    for (const Cell &cell: u.grid()->localActiveCells())
        u(cell) -= d(cell) * gradPCorr(cell);

    grid_->sendMessages(u); // gradPCorr may not be correct in buffer zones

    for (const Face &face: u.grid()->interiorFaces())
        u(face) -= d(face) * gradPCorr(face);

    for (const Face &face: u.grid()->boundaryFaces())
    {
        switch (u.boundaryType(face))
        {
            case VectorFiniteVolumeField::FIXED:
                break;

            case VectorFiniteVolumeField::SYMMETRY:
            {
                const Vector2D nWall = face.outwardNorm(face.lCell().centroid());
                u(face) = u(face.lCell()) - dot(u(face.lCell()), nWall) * nWall / nWall.magSqr();
            }
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                u(face) -= d(face) * gradPCorr(face);
                break;
        }
    }
}

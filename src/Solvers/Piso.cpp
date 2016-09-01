#include "Piso.h"
#include "Exception.h"
#include "CrankNicolson.h"
#include "AdamsBashforth.h"
#include "FaceInterpolation.h"
#include "GradientEvaluation.h"

Piso::Piso(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Solver(grid, input),
      u(addVectorField(input, "u")),
      h(addVectorField("h")),
      sg(addVectorField("sg")),
      gradP(addVectorField("gradP")),
      gradPCorr(addVectorField("gradPCorr")),
      p(addScalarField(input, "p")),
      pCorr(addScalarField("pCorr")),
      rho(addScalarField("rho")),
      mu(addScalarField("mu")),
      m(addScalarField("m")),
      d(addScalarField("d")),
      uEqn_(u, "momentum", SparseMatrix::IncompleteLUT),
      pCorrEqn_(pCorr, "pressure correction", SparseMatrix::IncompleteLUT)
{
    rho.fill(input.caseInput().get<Scalar>("Properties.rho", 1.));
    mu.fill(input.caseInput().get<Scalar>("Properties.mu", 1.));
    g_ = Vector2D(input.caseInput().get<std::string>("Properties.g", "(0,0)"));

    nInnerIterations_ = input.caseInput().get<size_t>("Solver.numInnerIterations");
    nPCorrections_ = input.caseInput().get<size_t>("Solver.numPressureCorrections");
    momentumOmega_ = input.caseInput().get<Scalar>("Solver.momentumRelaxation");
    pCorrOmega_ = input.caseInput().get<Scalar>("Solver.pressureCorrectionRelaxation");

    volumeIntegrators_ = VolumeIntegrator::initVolumeIntegrators(input, *this);

    forceIntegrators_ = ForceIntegrator::initForceIntegrators(input, p, rho, mu, u);

    pCorr.copyBoundaryTypes(p);

    uEqn_.matrix().setFill(1);
    pCorrEqn_.matrix().setFill(3);

    //- Perform a pseudo initialization
    u.savePreviousTimeStep(0., 1);
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

    for(const ForceIntegrator &fi: forceIntegrators_)
        fi.integrate();

    for(const VolumeIntegrator &vi: volumeIntegrators_)
        vi.integrate();

    printf("Max Co = %lf\n", courantNumber(timeStep));

    return 0.;
}

Scalar Piso::computeMaxTimeStep(Scalar maxCo) const
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

Scalar Piso::solveUEqn(Scalar timeStep)
{
    sg = fv::gravity(rho, g_);

    uEqn_ = (fv::ddt(rho, u, timeStep) + cn::div(rho*u, u) + ib_.eqns(u) == ab::laplacian(mu, u) - fv::source(gradP) + fv::source(sg));
    uEqn_.relax(momentumOmega_);

    Scalar error = uEqn_.solve();

    rhieChowInterpolation();

    return error;
}

Scalar Piso::solvePCorrEqn()
{
    for(const Cell& cell: m.grid.fluidCells())
    {
        m(cell) = 0.;

        for(const InteriorLink& nb: cell.neighbours())
            m(cell) += dot(u(nb.face()), nb.outwardNorm());

        for(const BoundaryLink& bd: cell.boundaries())
            m(cell) += dot(u(bd.face()), bd.outwardNorm());
    }

    pCorrEqn_ = (fv::laplacian(d, pCorr) + ib_.eqns(pCorr) == m);
    Scalar error = pCorrEqn_.solve();

    computeGradient(fv::FACE_TO_CELL, pCorr, gradPCorr);

    return error;
}

void Piso::rhieChowInterpolation()
{   
    VectorFiniteVolumeField& uStar = u.prevIter();
    VectorFiniteVolumeField& uPrev = u.prev(0);
    const Scalar dt = u.prevTimeStep(0);

    const auto diag = uEqn_.matrix().diagonal();

    d.fill(0.);
    for(const Cell& cell: d.grid.fluidCells())
        d(cell) = cell.volume()/diag[cell.globalIndex()];

    interpolateFaces(fv::INVERSE_VOLUME, d);

    for(const Face& face: u.grid.interiorFaces())
    {
        const Cell& cellP = face.lCell();
        const Cell& cellQ = face.rCell();

        const Scalar dP = d(cellP);
        const Scalar dQ = d(cellQ);
        const Scalar df = d(face);
        const Scalar rhoP = rho(cellP);
        const Scalar rhoQ = rho(cellQ);
        const Scalar rhof = rho(face);

        const Scalar g = cellQ.volume()/(cellP.volume() + cellQ.volume());

        if(!(cellP.isFluidCell() && cellQ.isFluidCell()))
        {
            u(face) = g*u(cellP) + (1. - g)*u(cellQ);
        }
        else
        {
            u(face) = g*u(cellP) + (1. - g)*u(cellQ)
                    + (1. - momentumOmega_)*(uStar(face) - (g*uStar(cellP) + (1. - g)*uStar(cellQ)))
                    + (rhof*df*uPrev(face) - (g*rhoP*dP*uPrev(cellP) + (1. - g)*rhoQ*dQ*uPrev(cellQ)))/dt //- This term is very important!
                    - df*gradP(face) + (g*dP*gradP(cellP) + (1. - g)*dQ*gradP(cellQ))
                    + df*sg(face) - (g*dP*sg(cellP) + (1. - g)*dQ*sg(cellQ));
        }
    }

    for(const Face& face: u.grid.boundaryFaces())
    {
        const Scalar df = d(face);

        switch(u.boundaryType(face))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            u(face) = u(face.lCell())
                    - df*gradP(face) + d(face.lCell())*gradP(face.lCell())
                    + df*sg(face) - d(face.lCell())*sg(face.lCell());
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

    interpolateFaces(fv::INVERSE_VOLUME, p);
    computeStaticPressure();
    computeGradient(fv::FACE_TO_CELL, p, gradP, true);
}

void Piso::correctVelocity()
{
    for(const Cell& cell: u.grid.fluidCells())
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

    for(const Cell& cell: m.grid.fluidCells())
    {
        m(cell) = 0.;

        for(const InteriorLink& nb: cell.neighbours())
            m(cell) += dot(u(nb.face()), nb.outwardNorm());

        for(const BoundaryLink& bd: cell.boundaries())
            m(cell) += dot(u(bd.face()), bd.outwardNorm());
    }
}

void Piso::computeStaticPressure()
{
    for(const Face& face: grid_.boundaryFaces())
    {
        if(p.boundaryType(face) == ScalarFiniteVolumeField::NORMAL_GRADIENT)
            p(face) = p(face.lCell()) + rho(face.lCell())*dot(g_, face.centroid() - face.lCell().centroid());
    }

    for(const ImmersedBoundaryObject& ibObj: ib_.ibObjs())
        for(const Cell &cell: ibObj.cells())
            for(const InteriorLink &nb: cell.neighbours())
            {
                if(grid_.fluidCells().isInGroup(nb.cell()) && ibObj.boundaryType(p.name) == ImmersedBoundaryObject::NORMAL_GRADIENT)
                {
                    p(nb.face()) = p(nb.cell()) + rho(nb.cell())*dot(g_, nb.face().centroid() - nb.cell().centroid());
                }
            }
}

Scalar Piso::courantNumber(Scalar timeStep)
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

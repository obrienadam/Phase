#include "ImmersedBoundary.h"
#include "Cicsam.h"
#include "GhostCellImmersedBoundary.h"
#include "CrankNicolson.h"
#include "AdamsBashforth.h"

ImmersedBoundary::ImmersedBoundary(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Multiphase(grid, input),
      csf_(input, gamma, u, scalarFields_, vectorFields_),
      cellStatus_(addScalarField("cell_status"))
{
    for(const auto& ibObject: input.boundaryInput().get_child("ImmersedBoundaries"))
    {
        printf("Initializing immersed boundary object \"%s\".\n", ibObject.first.c_str());

        ibObjs_.push_back(
                    ImmersedBoundaryObject(ibObject.first,
                                           grid,
                                           csf_,
                                           Vector2D(ibObject.second.get<std::string>("geometry.center")),
                                           ibObject.second.get<Scalar>("geometry.radius"))
                    );

        for(const auto& child: ibObject.second)
        {
            if(child.first == "geometry")
                continue;

            std::string type = child.second.get<std::string>("type");
            ImmersedBoundaryObject::BoundaryType boundaryType;
            Scalar boundaryRefValue = 0.;

            if(type == "fixed")
                boundaryType = ImmersedBoundaryObject::FIXED;
            else if(type == "normal_gradient")
                boundaryType = ImmersedBoundaryObject::NORMAL_GRADIENT;
            else if(type == "contact_angle")
                boundaryType = ImmersedBoundaryObject::CONTACT_ANGLE;
            else if(type == "partial_slip")
            {
                boundaryType = ImmersedBoundaryObject::PARTIAL_SLIP;
                boundaryRefValue = child.second.get<Scalar>("lambda");
            }
            else
                throw Exception("ImmersedBoundary", "ImmersedBoundary", "unrecognized boundary type \"" + type + "\".");

            printf("Setting boundary type \"%s\" for field \"%s\".\n", type.c_str(), child.first.c_str());

            ibObjs_.back().addBoundaryType(child.first, boundaryType);
            ibObjs_.back().addBoundaryRefValue(child.first, boundaryRefValue);

            if(child.first == "p")
                ibObjs_.back().addBoundaryType("pCorr", boundaryType);
        }

        ibObjs_.back().setInternalCells();
    }

    //- must reconstruct smoothing kernels
    csf_.constructSmoothingKernels();

    for(const Cell &cell: grid.inactiveCells())
    {
        u[cell.id()] = Vector2D(0., 0.);
        p[cell.id()] = 0.;
        gamma[cell.id()] = 0.;
    }

    // Only multiphase methods supported at the moment
    interfaceAdvectionMethod_ = CICSAM;
    setCellStatus();
}

Scalar ImmersedBoundary::solve(Scalar timeStep)
{
    computeRho();
    computeMu();

    u.savePreviousTimeStep(timeStep, 1);

    Scalar avgError = 0.;
    for(size_t innerIter = 0; innerIter < nInnerIterations_; ++innerIter)
    {
        u.savePreviousIteration();

        avgError += solveUEqn(timeStep);

        for(size_t pCorrIter = 0; pCorrIter < nPCorrections_; ++pCorrIter)
        {
            avgError += solvePCorrEqn();
            correctPressure();
            correctVelocity();
        }
    }

    solveGammaEqn(timeStep);

    printf("time step = %lf, max Co = %lf\n", timeStep, courantNumber(timeStep));

    return 0.;
}

Scalar ImmersedBoundary::solveUEqn(Scalar timeStep)
{
    sg = fv::gravity(rho, g_);
    ft = csf_.compute(ibObjs_);

    uEqn_ = (fv::ddt(rho, u, timeStep) + cn::div(rho*u, u) + gc::ib(ibObjs_, u)
             == ab::laplacian(mu, u) - fv::grad(p) + fv::source(ft) + fv::source(sg));

    uEqn_.relax(momentumOmega_);

    Scalar error = uEqn_.solve();

    rhieChowInterpolation();
    return error;
}

Scalar ImmersedBoundary::solvePCorrEqn()
{
    pCorrEqn_ = (fv::laplacian(d, pCorr) + gc::ib(ibObjs_, pCorr) == m);

    Scalar error = pCorrEqn_.solve();

    interpolateFaces(pCorr);

    return error;
}

Scalar ImmersedBoundary::solveGammaEqn(Scalar timeStep)
{
    gamma.savePreviousTimeStep(timeStep, 1);
    interpolateFaces(gamma);

    gammaEqn_ = (cicsam::div(u, gamma, timeStep, cicsam::HC) + gc::ib(ibObjs_, gamma) == 0.);

    Scalar error = gammaEqn_.solve();

    if(isnan(error))
        throw Exception("ImmersedBoundary", "solveGammaEqn", "a nan value was detected.");

    return error;
}

//- Protected

void ImmersedBoundary::setCellStatus()
{
    for(const Cell &cell: grid_.inactiveCells())
        cellStatus_[cell.id()] = SOLID;

    for(const Cell &cell: grid_.fluidCells())
        cellStatus_[cell.id()] = FLUID;

    for(const ImmersedBoundaryObject& ibObj: ibObjs_)
        for(const Cell &cell: ibObj.cells())
            cellStatus_[cell.id()] = IB;
}

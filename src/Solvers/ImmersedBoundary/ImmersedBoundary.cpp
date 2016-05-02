#include "ImmersedBoundary.h"
#include "Cicsam.h"
#include "GhostCellImmersedBoundary.h"

ImmersedBoundary::ImmersedBoundary(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Multiphase(grid, input),
      ibObj_(grid),
      cellStatus_(addScalarField("cell_status"))
{
    for(const auto& child: input.boundaryInput().get_child("ImmersedBoundaries"))
    {
        std::string name = child.first;

        const auto &ptree = child.second;

        ibObj_.init(Point2D(ptree.get<std::string>("geometry.center")),
                    ptree.get<Scalar>("geometry.radius"));

        ibObj_.addBoundaryType("u", ImmersedBoundaryObject::FIXED);
        ibObj_.addBoundaryType("p", ImmersedBoundaryObject::NORMAL_GRADIENT);
        ibObj_.addBoundaryType("pCorr", ImmersedBoundaryObject::NORMAL_GRADIENT);
        ibObj_.addBoundaryType("gamma", ImmersedBoundaryObject::NORMAL_GRADIENT);

        break; // for now only one ib object is allowed
    }

    ibObj_.constructStencils();
    setCellStatus();
}

Scalar ImmersedBoundary::solve(Scalar timeStep)
{
    computeRho();
    computeMu();

    u.save();

    Scalar avgError = 0.;
    for(size_t innerIter = 0; innerIter < nInnerIterations_; ++innerIter)
    {
        avgError += solveUEqn(timeStep);

        for(size_t pCorrIter = 0; pCorrIter < nPCorrections_; ++pCorrIter)
        {
            avgError += solvePCorrEqn();
            correctPressure();
            correctVelocity();
        }
    }

    printf("time step = %lf, max Co = %lf\n", timeStep, courantNumber(timeStep));

    return 0.;
}

Scalar ImmersedBoundary::solveUEqn(Scalar timeStep)
{
    computeInterfaceNormals();
    computeCurvature();

    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rho*u, u) + gc::ib(ibObj_, u)
             == fv::laplacian(mu, u) - fv::grad(p));

    uEqn_.relax(momentumOmega_);

    Scalar error = uEqn_.solve();

    rhieChowInterpolation();
    return error;
}

Scalar ImmersedBoundary::solvePCorrEqn()
{
    pCorrEqn_ = (fv::laplacian(d, pCorr) + gc::ib(ibObj_, pCorr) == m);

    Scalar error = pCorrEqn_.solve();

    interpolateFaces(pCorr);

    return error;
}

Scalar ImmersedBoundary::solveGammaEqn(Scalar timeStep)
{
    interpolateFaces(gamma);
    gamma.save();
    gammaEqn_ = (fv::ddt(gamma, timeStep) + cicsam::div(u, gamma, timeStep) + gc::ib(ibObj_, gamma) == 0.);

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

    for(const Cell &cell: grid_.cellGroup("ibCells"))
        cellStatus_[cell.id()] = IB;
}
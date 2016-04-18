#include "Multiphase.h"
#include "LineSegment2D.h"

Multiphase::Multiphase(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Piso(grid, input),
      gamma(addScalarField(input, "gamma")),
      kappa(addScalarField("kappa")),
      gammaEqn_(gamma, "gamma")
{
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1");
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2");
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1");
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2");
    sigma_ = input.caseInput().get<Scalar>("Properties.sigma");

    setInitialConditions(input);

    computeRho();
    computeMu();
}

Scalar Multiphase::solve(Scalar timeStep)
{
    computeRho();
    computeMu();

    Piso::solve(timeStep);

    solveGammaEqn(timeStep);
}

//- Private methods

void Multiphase::computeRho()
{
    for(const Cell& cell: rho.grid.cells())
    {
        size_t id = cell.id();
        rho[id] = rho1_*(1. - gamma[id]) + rho2_*gamma[id];
    }

    interpolateFaces(rho);
}

void Multiphase::computeMu()
{
    for(const Cell& cell: mu.grid.cells())
    {
        size_t id = cell.id();
        mu[id] = mu1_*(1. - gamma[id]) + mu2_*gamma[id];
    }

    interpolateFaces(mu);
}

Scalar Multiphase::solveGammaEqn(Scalar timeStep)
{
    interpolateFaces(gamma);
    gamma.save();
    gammaEqn_ = (fv::ddt(gamma, timeStep) + cicsam::div(u, gamma, timeStep) == 0.);

    Scalar error = gammaEqn_.solve();

    if(isnan(error))
        throw Exception("Multiphase", "solveGammaEqn", "a nan value was detected.");

    return error;
}

//- External functions

namespace hc
{

Scalar betaFace(Scalar gammaD, Scalar gammaA, Scalar gammaU, Scalar coD)
{
    const Scalar gammaTilde = (gammaD - gammaU)/(gammaA - gammaU);
    const Scalar gammaTildeFace = gammaTilde >= 0 && gammaTilde <= 1 ? std::min(1., gammaTilde/coD) : gammaTilde;
    const Scalar betaFace = (gammaTildeFace - gammaTilde)/(1. - gammaTilde);

    return betaFace >= 0. && betaFace <= 1. ? betaFace : 0.5;
}

}

namespace sc
{

Scalar betaFace(Scalar gammaD, Scalar gammaA, Scalar gammaU, Scalar coD)
{
    const Scalar gammaTilde = (gammaD - gammaU)/(gammaA - gammaU);
    const Scalar gammaTildeFace = gammaTilde >= 0 && gammaTilde <= 1 ? std::min(((1. - coD) + (1. + coD)*gammaTilde)/2., std::min(1., gammaTilde/coD)) : gammaTilde;
    const Scalar betaFace = (gammaTildeFace - gammaTilde)/(1. - gammaTilde);

    return betaFace >= 0. && betaFace <= 1. ? betaFace : 0.5;
}

}

namespace uq
{

Scalar betaFace(Scalar gammaD, Scalar gammaA, Scalar gammaU, Scalar coD)
{
    const Scalar gammaTilde = (gammaD - gammaU)/(gammaA - gammaU);
    const Scalar gammaTildeFace = gammaTilde >= 0 && gammaTilde <= 1 ? std::min((8.*coD*gammaTilde + (1. - coD)*(6.*gammaTilde + 3.))/8., std::min(1., gammaTilde/coD)) : gammaTilde;
    const Scalar betaFace = (gammaTildeFace - gammaTilde)/(1. - gammaTilde);

    return betaFace >= 0. && betaFace <= 1. ? betaFace : 0.5;
}

}

namespace cicsam
{

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField &u, ScalarFiniteVolumeField &field, Scalar timeStep)
{
    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(field);
    VectorFiniteVolumeField gradField = grad(field);

    entries.reserve(5*field.grid.nActiveCells());

    for(const Cell &cell: field.grid.cells())
    {
        Scalar centralCoeff = 0.;
        size_t row = cell.globalIndex();

        for(const InteriorLink &nb: cell.neighbours())
        {
            //- Determine the donor and acceptor cell

            const Scalar flux = dot(u.faces()[nb.face().id()], nb.outwardNorm());
            const Cell &donor = flux > 0. ? cell : nb.cell();
            const Cell &acceptor = flux > 0. ? nb.cell() : cell;

            const Scalar gammaD = field[donor.id()];
            const Scalar gammaA = field[acceptor.id()];
            Scalar gammaU = gammaA + 2*dot(donor.centroid() - acceptor.centroid(), gradField[donor.id()]);

            if(gammaU < 0)
                gammaU = 0;
            else if(gammaU > 1)
                gammaU = 1;


            const Scalar coD = u.faces()[nb.face().id()].mag()*timeStep/nb.rCellVec().mag();

            Scalar betaFace = uq::betaFace(gammaD, gammaA, gammaU, coD);

            size_t col = nb.cell().globalIndex();
            Scalar coeff;

            if(&cell == &donor)
            {
                coeff = betaFace*flux;
                centralCoeff += (1. - betaFace)*flux;
            }
            else
            {
                coeff = (1. - betaFace)*flux;
                centralCoeff += betaFace*flux;
            }

            entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, col, coeff));
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Scalar flux = dot(u.faces()[bd.face().id()], bd.outwardNorm());

            switch(field.boundaryType(bd.face().id()))
            {
            case ScalarFiniteVolumeField::FIXED:
                eqn.boundaries()(row) -= flux*field.faces()[bd.face().id()];
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                centralCoeff += flux;
                break;

            default:
                throw Exception("hc", "div", "unrecognized or unspecified boundary type.");
            }
        }

        entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, row, centralCoeff));
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

} // end namespace cicsam

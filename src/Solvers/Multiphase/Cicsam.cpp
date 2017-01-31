#include "Cicsam.h"

namespace cicsam
{

Scalar hc(Scalar gammaTilde, Scalar coD)
{
    return gammaTilde >= 0 && gammaTilde <= 1 ? std::min(1., gammaTilde/coD) : gammaTilde;
}

Scalar uq(Scalar gammaTilde, Scalar coD)
{
    return gammaTilde >= 0 && gammaTilde <= 1 ? std::min((8.*coD*gammaTilde + (1. - coD)*(6.*gammaTilde + 3.))/8., hc(gammaTilde, coD)) : gammaTilde;
}

Equation<Scalar> cn(const VectorFiniteVolumeField &u,
                    const VectorFiniteVolumeField& gradGamma,
                    const VectorFiniteVolumeField &m,
                    ScalarFiniteVolumeField &gamma,
                    Scalar timeStep)
{
    Equation<Scalar> eqn(gamma);
    const Scalar k = 1; //- 0 For a full UQ scheme, 2 for max HC

    for(const Cell &cell: gamma.grid.cellZone("fluid"))
    {
        Scalar centralCoeff = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            //- Determine the donor and acceptor cell

            Scalar flux = dot(u(nb.face()), nb.outwardNorm());
            const Cell &donor = flux > 0. ? cell : nb.cell();
            const Cell &acceptor = flux > 0. ? nb.cell() : cell;

            const Scalar gammaD = gamma(donor);
            const Scalar gammaA = gamma(acceptor);
            const Scalar gammaU = std::max(std::min(gammaA + 2*dot(donor.centroid() - acceptor.centroid(), gradGamma(donor)), 1.), 0.);

            //const Scalar coD = sqrt(u(nb.face()).magSqr()*timeStep*timeStep/nb.rCellVec().magSqr());
            const Scalar coD = fabs(dot(u(nb.face()), nb.outwardNorm())/dot(nb.rCellVec(), nb.outwardNorm())*timeStep);

            const Scalar gammaTilde = (gammaD - gammaU)/(gammaA - gammaU);

            const Scalar thetaF = acos(fabs(dot(m(nb.face()), nb.rCellVec().unitVec())));
            const Scalar psiF = std::min(k*(cos(2*thetaF) + 1.)/2., 1.);

            const Scalar gammaTildeF = psiF*hc(gammaTilde, coD) + (1. - psiF)*uq(gammaTilde, coD);

            Scalar betaFace = (gammaTildeF - gammaTilde)/(1. - gammaTilde);

            if(std::isnan(betaFace)) // Central difference if non-valid beta value
                betaFace = &cell == &donor ? 0. : 1.;

            betaFace = std::max(std::min(betaFace, 1.), 0.);

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

            eqn.add(cell, nb.cell(), coeff/2.);
            eqn.addBoundary(cell, -coeff*gamma.prev()(nb.cell())/2.);
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Scalar flux = dot(u(bd.face()), bd.outwardNorm());

            switch(gamma.boundaryType(bd.face()))
            {
            case ScalarFiniteVolumeField::FIXED:
                eqn.addBoundary(cell, -flux*gamma(bd.face())/2.);
                eqn.addBoundary(cell, -flux*gamma.prev(0)(bd.face())/2.);
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT: case ScalarFiniteVolumeField::SYMMETRY:
                centralCoeff += flux;
                break;

            default:
                throw Exception("cicsam", "cn", "unrecognized or unspecified boundary type.");
            }
        }

        eqn.add(cell, cell, centralCoeff/2.);
        eqn.addBoundary(cell, -centralCoeff*gamma.prev()(cell)/2.);
    }

    return eqn;
}

} // end namespace cicsam

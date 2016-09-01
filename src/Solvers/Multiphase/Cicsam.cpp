#include "Cicsam.h"

namespace cicsam
{

Scalar hc(Scalar gammaTilde, Scalar coD) // Maintains very sharp interface but funny things sometimes happen
{
    return gammaTilde >= 0 && gammaTilde <= 1 ? std::min(1., gammaTilde/coD) : gammaTilde;
}


Scalar uq(Scalar gammaTilde, Scalar coD) // Fairly sharp interface but fewer funny things
{
    return gammaTilde >= 0 && gammaTilde <= 1 ? std::min((8.*coD*gammaTilde + (1. - coD)*(6.*gammaTilde + 3.))/8., hc(gammaTilde, coD)) : gammaTilde;
}

Equation<ScalarFiniteVolumeField> cn(const VectorFiniteVolumeField &u,
                                     const VectorFiniteVolumeField& gradGamma,
                                     const VectorFiniteVolumeField &m,
                                     ScalarFiniteVolumeField &gamma,
                                     Scalar timeStep)
{
    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(gamma);
    const Scalar k = 1.5; //- 0 For a full UQ scheme, 2 for max HC

    entries.reserve(5*gamma.grid.nActiveCells());

    for(const Cell &cell: gamma.grid.fluidCells())
    {
        Scalar centralCoeff = 0.;
        const Index row = cell.globalIndex();

        for(const InteriorLink &nb: cell.neighbours())
        {
            //- Determine the donor and acceptor cell

            const Scalar flux = dot(u(nb.face()), nb.outwardNorm());
            const Cell &donor = flux > 0. ? cell : nb.cell();
            const Cell &acceptor = flux > 0. ? nb.cell() : cell;

            const Scalar gammaD = gamma(donor);
            const Scalar gammaA = gamma(acceptor);
            const Scalar gammaU = std::max(std::min(gammaA + 2*dot(donor.centroid() - acceptor.centroid(), gradGamma(donor)), 1.), 0.);

            const Scalar coD = sqrt(u(nb.face()).magSqr()*timeStep*timeStep/nb.rCellVec().magSqr());

            const Scalar gammaTilde = (gammaD - gammaU)/(gammaA - gammaU);

            const Scalar thetaF = acos(fabs(dot(m(cell), nb.rCellVec().unitVec())));
            const Scalar psiF = std::min(k*(cos(2*thetaF) + 1.)/2., 1.);

            const Scalar gammaTildeF = psiF*hc(gammaTilde, coD) + (1. - psiF)*uq(gammaTilde, coD);

            Scalar betaFace = (gammaTildeF - gammaTilde)/(1. - gammaTilde);

            if(isnan(betaFace)) // Upwind if a valid value is not found
                betaFace = &cell == &donor ? 0. : 1.;

            betaFace = std::max(std::min(betaFace, 1.), 0.);

            const size_t col = nb.cell().globalIndex();
            Scalar coeff;
            if(&cell == &donor)
            {
                coeff = betaFace*flux/2.;
                centralCoeff += (1. - betaFace)*flux/2.;
            }
            else
            {
                coeff = (1. - betaFace)*flux/2.;
                centralCoeff += betaFace*flux/2.;
            }

            entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, col, coeff));
            eqn.boundaries()(row) -= coeff*gamma.prev()(nb.cell());
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Scalar flux = dot(u(bd.face()), bd.outwardNorm());

            switch(gamma.boundaryType(bd.face()))
            {
            case ScalarFiniteVolumeField::FIXED:
                eqn.boundaries()(row) -= flux*gamma(bd.face());
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                centralCoeff += flux/2.;
                eqn.boundaries()(row) -= flux/2.*gamma.prev()(cell);
                break;

            case ScalarFiniteVolumeField::SYMMETRY:
                break;

            default:
                throw Exception("cicsam", "cn", "unrecognized or unspecified boundary type.");
            }
        }

        entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, row, centralCoeff));
        eqn.boundaries()(row) -= centralCoeff*gamma.prev()(cell);
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

} // end namespace cicsam

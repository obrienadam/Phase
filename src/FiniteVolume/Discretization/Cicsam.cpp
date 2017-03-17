#include "Cicsam.h"

namespace cicsam
{

Scalar hc(Scalar gammaTildeD, Scalar coD)
{
    return gammaTildeD >= 0 && gammaTildeD <= 1 ? std::min(1., gammaTildeD / coD) : gammaTildeD;
}

Scalar uq(Scalar gammaTildeD, Scalar coD)
{
    return gammaTildeD >= 0 && gammaTildeD <= 1 ? std::min(
                                                      (8. * coD * gammaTildeD + (1. - coD) * (6. * gammaTildeD + 3.)) / 8., hc(gammaTildeD, coD)) : gammaTildeD;
}

Equation<Scalar> cn(const VectorFiniteVolumeField &u,
                    const VectorFiniteVolumeField &gradGamma,
                    const VectorFiniteVolumeField &m,
                    ScalarFiniteVolumeField &gamma,
                    Scalar timeStep, Scalar alpha, Scalar k)
{
    Equation<Scalar> eqn(gamma);
    const VectorFiniteVolumeField& u0 = u.prev(0);
    const CellZone& fluid = gamma.grid.cellZone("fluid");

    for (const Cell &cell: fluid)
    {
        Scalar centralCoeff = 0., centralCoeff0 = 0.;

        for (const InteriorLink &nb: cell.neighbours())
        {
            //- Determine the donor and acceptor cell

            Scalar flux = dot(u(nb.face()), nb.outwardNorm());
            Scalar flux0 = dot(u0(nb.face()), nb.outwardNorm());

            const Cell &donor = flux > 0. ? cell : nb.cell();
            const Cell &acceptor = flux > 0. ? nb.cell() : cell;

            const Scalar gammaD = gamma(donor);
            const Scalar gammaA = gamma(acceptor);
            const Scalar gammaU = std::max(
                        std::min(gammaA + 2 * dot(donor.centroid() - acceptor.centroid(), gradGamma(donor)), 1.), 0.);

            //Scalar gammaU = gamma(fluid.cellNearestNeighbours(2*donor.centroid() - acceptor.centroid(), 1)[0]);

            //const Scalar coD = sqrt(u(nb.face()).magSqr()*timeStep*timeStep/nb.rCellVec().magSqr());
            const Scalar coD = fabs(
                        dot(u(nb.face()), nb.outwardNorm()) / dot(nb.rCellVec(), nb.outwardNorm()) * timeStep);

            const Scalar gammaTilde = (gammaD - gammaU) / (gammaA - gammaU);

            Scalar thetaF = acos(fabs(dot(m(nb.face()), nb.rCellVec().unitVec())));
            Scalar psiF = std::min(k * (cos(2 * thetaF) + 1.) / 2., 1.);
            Scalar gammaTildeF = psiF * hc(gammaTilde, coD) + (1. - psiF) * uq(gammaTilde, coD);
            Scalar betaFace = (gammaTildeF - gammaTilde) / (1. - gammaTilde);

            if (std::isnan(betaFace))
                betaFace = &cell == &donor ? 0. : 1.;

            betaFace = std::max(std::min(betaFace, 1.), 0.);

            Scalar coeff, coeff0;
            if (cell.id() == donor.id())
            {
                coeff = betaFace * flux;
                centralCoeff += (1. - betaFace) * flux;

                coeff0 = betaFace * flux0;
                centralCoeff0 += (1. - betaFace) * flux0;
            }
            else
            {
                coeff = (1. - betaFace) * flux;
                centralCoeff += betaFace * flux;

                coeff0 = (1. - betaFace) * flux0;
                centralCoeff0 += betaFace * flux0;
            }

            eqn.add(cell, nb.cell(), coeff * alpha);
            eqn.addBoundary(cell, -coeff0 * gamma.prev()(nb.cell()) * (1. - alpha));
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Scalar flux = dot(u(bd.face()), bd.outwardNorm());
            Scalar flux0 = dot(u0(bd.face()), bd.outwardNorm());

            switch (gamma.boundaryType(bd.face()))
            {
            case ScalarFiniteVolumeField::FIXED:
                eqn.addBoundary(cell, -flux * gamma(bd.face()) * alpha);
                eqn.addBoundary(cell, -flux0 * gamma.prev(0)(bd.face()) * (1. - alpha));
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
            case ScalarFiniteVolumeField::SYMMETRY:
                centralCoeff += flux;
                centralCoeff0 += flux0;
                break;

            default:
                throw Exception("cicsam", "cn", "unrecognized or unspecified boundary type.");
            }
        }

        eqn.add(cell, cell, centralCoeff * alpha);
        eqn.addBoundary(cell, -centralCoeff0 * gamma.prev()(cell) * (1. - alpha));
    }

    return eqn;
}

void interpolateFaces(const VectorFiniteVolumeField &u,
                      const VectorFiniteVolumeField &gradGamma,
                      const VectorFiniteVolumeField &m,
                      ScalarFiniteVolumeField &gamma,
                      Scalar timeStep, Scalar k)
{
    for (const Face &face: gamma.grid.interiorFaces())
    {
        //- Determine the donor and acceptor cell
        Vector2D sf = face.outwardNorm(face.lCell().centroid());
        Scalar flux = dot(u(face), sf);

        const Cell &donor = flux > 0. ? face.lCell() : face.rCell();
        const Cell &acceptor = flux > 0. ? face.rCell() : face.lCell();

        const Scalar gammaD = gamma(donor);
        const Scalar gammaA = gamma(acceptor);
        const Scalar gammaU = std::max(
                    std::min(gammaA + 2 * dot(donor.centroid() - acceptor.centroid(), gradGamma(donor)), 1.), 0.);

        const Scalar coD = fabs(
                    dot(u(face), sf) / dot(acceptor.centroid() - donor.centroid(), sf) * timeStep);

        const Scalar gammaTilde = (gammaD - gammaU) / (gammaA - gammaU);

        Scalar thetaF = acos(fabs(dot(m(face), (donor.centroid() - acceptor.centroid()).unitVec())));
        Scalar psiF = std::min(k * (cos(2 * thetaF) + 1.) / 2., 1.);
        Scalar gammaTildeF = psiF * hc(gammaTilde, coD) + (1. - psiF) * uq(gammaTilde, coD);
        Scalar betaFace = (gammaTildeF - gammaTilde) / (1. - gammaTilde);

        if (std::isnan(betaFace))
            betaFace = face.lCell().id() == donor.id() ? 0. : 1.;

        betaFace = std::max(std::min(betaFace, 1.), 0.);

        if (face.lCell().id() == donor.id())
            gamma(face) = (1. - betaFace)*gamma(face.lCell()) + betaFace*gamma(face.rCell());
        else
            gamma(face) = betaFace*gamma(face.lCell()) + (1. - betaFace)*gamma(face.rCell());
    }

    for (const Face& face: gamma.grid.boundaryFaces())
    {
        switch (gamma.boundaryType(face))
        {
        case ScalarFiniteVolumeField::FIXED:
            break;

        case ScalarFiniteVolumeField::NORMAL_GRADIENT:
        case ScalarFiniteVolumeField::SYMMETRY:
            gamma(face) = gamma(face.lCell());
            break;

        default:
            throw Exception("cicsam", "interpolateFaces", "unrecognized or unspecified boundary type.");
        }
    }
}

} // end namespace cicsam

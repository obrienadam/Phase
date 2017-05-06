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

void interpolateFaces(const VectorFiniteVolumeField &u,
                      const VectorFiniteVolumeField &gradGamma,
                      ScalarFiniteVolumeField &gamma,
                      Scalar timeStep,
                      Scalar k)
{
    for(const Face& face: gamma.grid.interiorFaces())
    {
        Vector2D sf = face.outwardNorm(face.lCell().centroid());
        Scalar flux = dot(u(face), sf);

        const Cell& donor = flux >= 0. ? face.lCell() : face.rCell();
        const Cell& acceptor = flux >= 0. ? face.rCell() : face.lCell();
        Vector2D rc = acceptor.centroid() - donor.centroid();

        Scalar gammaD = gamma(donor);
        Scalar gammaA = gamma(acceptor);
        Scalar gammaU = std::max(std::min(gammaA - 2.*dot(rc, gradGamma(donor)), 1.), 0.);
        Scalar coD = fabs(dot(u(face), sf) / dot(acceptor.centroid() - donor.centroid(), sf) * timeStep);
        Scalar gammaTilde = (gammaD - gammaU) / (gammaA - gammaU);

        Scalar thetaF = acos(dot(gradGamma(donor), rc)*dot(gradGamma(donor), rc)/(gradGamma(donor).magSqr()*rc.magSqr()));
        Scalar psiF = std::min(k * (cos(2 * thetaF) + 1.) / 2., 1.);
        Scalar gammaTildeF = psiF * hc(gammaTilde, coD) + (1. - psiF) * uq(gammaTilde, coD);
        Scalar betaFace = (gammaTildeF - gammaTilde) / (1. - gammaTilde);

        if (std::isnan(betaFace))
            betaFace = 0.; // Default to upwind, the most stable option

        //- Make sure the stencil is bounded
        betaFace = std::max(std::min(betaFace, 1.), 0.);
        gamma(face) = (1. - betaFace)*gamma(donor) + betaFace*gamma(acceptor);
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

Equation<Scalar> div(const VectorFiniteVolumeField &u,
                     ScalarFiniteVolumeField &gamma)
{
    Equation<Scalar> eqn(gamma);

    for (const Cell &cell: gamma.grid.cellZone("fluid"))
    {
        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar flux = dot(u(nb.face()), nb.outwardNorm());
            eqn.addSource(cell, flux*gamma(nb.face()));
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Scalar flux = dot(u(bd.face()), bd.outwardNorm());
            switch (gamma.boundaryType(bd.face()))
            {
            case ScalarFiniteVolumeField::FIXED:
                eqn.addSource(cell, flux * gamma(bd.face()));
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                eqn.add(cell, cell, flux);
                break;

            case ScalarFiniteVolumeField::SYMMETRY:
                break;

            default:
                throw Exception("cicsam", "div", "unrecognized or unspecified boundary type.");
            }
        }
    }

    return eqn;
}

Equation<Scalar> div(const VectorFiniteVolumeField &u,
                     const VectorFiniteVolumeField &gradGamma,
                     ScalarFiniteVolumeField &gamma,
                     Scalar timeStep,
                     Scalar k)
{
    Equation<Scalar> eqn(gamma);

    for(const Cell& cell: gamma.grid.cellZone("fluid"))
    {
        for(const InteriorLink& nb: cell.neighbours())
        {
            Scalar flux = dot(u(nb.face()), nb.outwardNorm());

            if(flux < 0.) //- Not a donor cell for this face
                continue;

            const Cell& donor = cell;
            const Cell& acceptor = nb.cell();
            Vector2D rc = acceptor.centroid() - donor.centroid();

            Scalar gammaD = gamma(donor);
            Scalar gammaA = gamma(acceptor);
            Scalar gammaU = std::max(std::min(gammaA - 2.*dot(rc, gradGamma(donor)), 1.), 0.);
            Scalar coD = fabs(flux / dot(acceptor.centroid() - donor.centroid(), nb.outwardNorm()) * timeStep);
            Scalar gammaTilde = (gammaD - gammaU) / (gammaA - gammaU);

            Scalar thetaF = acos(dot(gradGamma(donor), rc)*dot(gradGamma(donor), rc)/(gradGamma(donor).magSqr()*rc.magSqr()));
            Scalar psiF = std::min(k * (cos(2 * thetaF) + 1.) / 2., 1.);
            Scalar gammaTildeF = psiF * hc(gammaTilde, coD) + (1. - psiF) * uq(gammaTilde, coD);
            Scalar betaFace = (gammaTildeF - gammaTilde) / (1. - gammaTilde);

            if (std::isnan(betaFace))
                betaFace = 0.; // Default to upwind, the most stable option

            //- Make sure the stencil is bounded
            betaFace = std::max(std::min(betaFace, 1.), 0.);

            eqn.add(donor, donor, flux*(1. - betaFace));
            eqn.add(donor, acceptor, flux*betaFace);

            eqn.add(acceptor, acceptor, -flux*betaFace);
            eqn.add(acceptor, donor, -flux*(1. - betaFace));
        }

        for (const BoundaryLink& bd: cell.boundaries())
        {
            Scalar flux = dot(u(bd.face()), bd.outwardNorm());

            switch (gamma.boundaryType(bd.face()))
            {
            case ScalarFiniteVolumeField::FIXED:
                eqn.addSource(cell, flux * gamma(bd.face()));
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                eqn.add(cell, cell, flux);
                break;

            case ScalarFiniteVolumeField::SYMMETRY:
                break;

            default:
                throw Exception("cicsam", "div", "unrecognized or unspecified boundary type.");
            }
        }
    }

    return eqn;
}

} // end namespace cicsam

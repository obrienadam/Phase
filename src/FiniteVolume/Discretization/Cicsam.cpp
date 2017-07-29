#include "Cicsam.h"
#include "Algorithm.h"

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

    ScalarFiniteVolumeField computeBeta(const VectorFiniteVolumeField &u,
                                        const VectorFiniteVolumeField &gradGamma,
                                        const ScalarFiniteVolumeField &gamma,
                                        Scalar timeStep,
                                        Scalar k)
    {
        ScalarFiniteVolumeField beta(gamma.gridPtr(), "beta");

        for(const Face& face: gamma.grid().interiorFaces())
        {
            Vector2D sf = face.outwardNorm(face.lCell().centroid());
            Scalar flux = dot(u(face), sf);

            const Cell& donor = flux >= 0. ? face.lCell() : face.rCell();
            const Cell& acceptor = flux >= 0. ? face.rCell() : face.lCell();
            Vector2D rc = acceptor.centroid() - donor.centroid();

            Scalar gammaD = clamp(gamma(donor), 0., 1.);
            Scalar gammaA = clamp(gamma(acceptor), 0., 1.);
            Scalar gammaU = clamp(gammaA - 2.*dot(rc, gradGamma(donor)), 0., 1.);
            Scalar gammaTilde = (gammaD - gammaU) / (gammaA - gammaU);
            Scalar coD = std::abs(dot(u(face), sf) / dot(rc, sf) * timeStep);
            Scalar thetaF = acos(fabs(dot((gradGamma(donor)/2 + gradGamma(acceptor)/2).unitVec(), rc.unitVec())));
            Scalar psiF = std::min(k * (cos(2 * thetaF) + 1.) / 2., 1.);
            Scalar gammaTildeF = psiF * hc(gammaTilde, coD) + (1. - psiF) * uq(gammaTilde, coD);
            Scalar betaFace = (gammaTildeF - gammaTilde) / (1. - gammaTilde);

            if (std::isnan(betaFace))
                betaFace = 0.5; // Default to a centered scheme

            //- Make sure the stencil is bounded
            beta(face) = clamp(betaFace, 0., 1.);
        }

        return beta;
    }

    Equation<Scalar> div(const VectorFiniteVolumeField &u,
                         const ScalarFiniteVolumeField &beta,
                         ScalarFiniteVolumeField &gamma,
                         Scalar theta)
    {
        Equation<Scalar> eqn(gamma);

        for (const Cell &cell: gamma.grid().cellZone("fluid"))
        {
            for (const InteriorLink &nb: cell.neighbours())
            {
                Scalar flux = dot(u(nb.face()), nb.outwardNorm());
                Scalar alpha = flux > 0 ? 1. - beta(nb.face()): beta(nb.face());

                eqn.add(cell, cell, theta*alpha*flux);
                eqn.add(cell, nb.cell(), theta*(1. - alpha)*flux);
                eqn.addSource(cell, (1. - theta)*flux*gamma(nb.face()));
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
                        eqn.add(cell, cell, theta*flux);
                        eqn.addSource(cell, (1. - theta)*flux*gamma(bd.face()));
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

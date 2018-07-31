#include "Cicsam.h"

#include "Math/Algorithm.h"

std::vector<Scalar> cicsam::beta(const VectorFiniteVolumeField &u,
                                 const VectorFiniteVolumeField &gradGamma,
                                 const ScalarFiniteVolumeField &gamma,
                                 Scalar timeStep,
                                 Scalar k)
{
    auto hc = [](Scalar gammaDTilde, Scalar coD) {
        return gammaDTilde >= 0 && gammaDTilde <= 1 ? std::min(1., gammaDTilde / coD) : gammaDTilde;
    };

    auto uq = [&hc](Scalar gammaDTilde, Scalar coD) {
        return gammaDTilde >= 0 && gammaDTilde <= 1 ?
                    std::min((8. * coD * gammaDTilde + (1. - coD) * (6. * gammaDTilde + 3.)) / 8., hc(gammaDTilde, coD)) :
                    gammaDTilde;
    };

    std::vector<Scalar> beta(gamma.grid()->faces().size(), 0.);

    for (const Face &face: gamma.grid()->interiorFaces())
    {
        Vector2D sf = face.outwardNorm(face.lCell().centroid());
        Scalar flux = dot(u(face), sf);
        const Cell &donor = flux >= 0. ? face.lCell() : face.rCell();
        const Cell &acceptor = flux >= 0. ? face.rCell() : face.lCell();
        //const Cell &upwind = gamma.grid()->globalCells().nearestItem(2. * donor.centroid() - acceptor.centroid());

        Vector2D rc = acceptor.centroid() - donor.centroid();

        Scalar gammaD = clamp(gamma(donor), 0., 1.);
        Scalar gammaA = clamp(gamma(acceptor), 0., 1.);
        Scalar gammaU = clamp(gammaA - 2. * dot(rc, gradGamma(donor)), 0., 1.);
        // Scalar gammaU = clamp(gamma(upwind), 0., 1.);

        Scalar gammaDTilde = (gammaD - gammaU) / (gammaA - gammaU);

        Scalar coD = 0.; //- Cell courant number
        for (const InteriorLink &nb: donor.neighbours())
            coD += std::max(dot(u(nb.face()), nb.outwardNorm()) / donor.volume() * timeStep, 0.);

        for (const BoundaryLink &bd: donor.boundaries())
            coD += std::max(dot(u(bd.face()), bd.outwardNorm()) / donor.volume() * timeStep, 0.);

        Scalar thetaF = std::acos(std::abs(dot(gradGamma(donor).unitVec(), rc.unitVec())));
        Scalar psiF = std::min(k * (std::cos(2 * thetaF) + 1.) / 2., 1.);
        Scalar gammaFTilde = psiF * hc(gammaDTilde, coD) + (1. - psiF) * uq(gammaDTilde, coD);
        Scalar betaFace = (gammaFTilde - gammaDTilde) / (1. - gammaDTilde);

        //- If stencil cannot be computed, default to upwind
        beta[face.id()] = std::isnan(betaFace) ? 0. : clamp(betaFace, 0., 1.);
    }

    return beta;
}

void cicsam::computeMomentumFlux(Scalar rho1,
                                 Scalar rho2,
                                 const VectorFiniteVolumeField &u,
                                 const ScalarFiniteVolumeField &gamma,
                                 const std::vector<Scalar> &beta,
                                 Scalar timeStep,
                                 VectorFiniteVolumeField &rhoU)
{
    rhoU.savePreviousTimeStep(timeStep, 1);

    rhoU.oldField(0).computeInteriorFaces([rho1, rho2, &u, &gamma, &beta](const Face &f) {
        Scalar flux = dot(u(f), f.outwardNorm());
        const Cell &d = flux > 0. ? f.lCell() : f.rCell();
        const Cell &a = flux > 0. ? f.rCell() : f.lCell();
        Scalar g = (1. - beta[f.id()]) * gamma.oldField(0)(d) + beta[f.id()] * gamma.oldField(0)(a);
        return ((1. - g) * rho1 + g * rho2) * u(f);
    });

    rhoU.oldField(0).computeBoundaryFaces([rho1, rho2, &u, &gamma](const Face &f) {
        Scalar g = gamma.oldField(0)(f);
        return ((1. - g) * rho1 + g * rho2) * u(f);
    });

    rhoU.computeInteriorFaces([rho1, rho2, &u, &gamma, &beta](const Face &f) {
        Scalar flux = dot(u(f), f.outwardNorm());
        const Cell &d = flux > 0. ? f.lCell() : f.rCell();
        const Cell &a = flux > 0. ? f.rCell() : f.lCell();
        Scalar g = (1. - beta[f.id()]) * gamma(d) + beta[f.id()] * gamma(a);
        return ((1. - g) * rho1 + g * rho2) * u(f);
    });

    rhoU.computeBoundaryFaces([rho1, rho2, &u, &gamma](const Face &f) {
        Scalar g = gamma(f);
        return ((1. - g) * rho1 + g * rho2) * u(f);
    });
}

FiniteVolumeEquation<Scalar> cicsam::div(const VectorFiniteVolumeField &u,
                                         const std::vector<Scalar> &beta,
                                         ScalarFiniteVolumeField &gamma,
                                         Scalar theta,
                                         const CellGroup &cells)
{
    FiniteVolumeEquation<Scalar> eqn(gamma);

    for (const Cell &cell: cells)
    {
        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar flux = dot(u(nb.face()), nb.outwardNorm());

            const Cell &donor = flux > 0. ? cell : nb.cell();
            const Cell &acceptor = flux > 0. ? nb.cell() : cell;

            Scalar b = beta[nb.face().id()];

            eqn.add(cell, donor, theta * (1. - b) * flux);
            eqn.add(cell, acceptor, theta * b * flux);

            Scalar gammaF = (1. - b) * gamma(donor) + b * gamma(acceptor);
            eqn.addSource(cell, (1. - theta) * flux * gammaF);
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
                eqn.add(cell, cell, theta * flux);
                eqn.addSource(cell, (1. - theta) * flux * gamma(bd.face()));
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

FiniteVolumeEquation<Scalar> cicsam::div(const VectorFiniteVolumeField &u,
                                         const std::vector<Scalar> &beta,
                                         ScalarFiniteVolumeField &gamma,
                                         Scalar theta)
{
    return div(u, beta, gamma, theta, gamma.cells());
}

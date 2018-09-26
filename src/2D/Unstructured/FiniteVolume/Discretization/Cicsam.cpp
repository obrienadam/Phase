#include "Cicsam.h"

#include "Math/Algorithm.h"

Scalar cicsam::hc(Scalar gammaDTilde, Scalar coD)
{
    return gammaDTilde >= 0 && gammaDTilde <= 1 ? std::min(1., gammaDTilde / coD) : gammaDTilde;
}

Scalar cicsam::uq(Scalar gammaDTilde, Scalar coD)
{
    return gammaDTilde >= 0 && gammaDTilde <= 1 ?
                std::min((8. * coD * gammaDTilde + (1. - coD) * (6. * gammaDTilde + 3.)) / 8., hc(gammaDTilde, coD)) :
                gammaDTilde;
}

std::vector<Scalar> cicsam::faceInterpolationWeights(const VectorFiniteVolumeField &u,
                                                     const ScalarFiniteVolumeField &gamma,
                                                     const VectorFiniteVolumeField &gradGamma,
                                                     Scalar timeStep)
{
    const Scalar k = 1;

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
                                 const std::vector<Scalar> &faceInterpolationWeights,
                                 VectorFiniteVolumeField &rhoU)
{
    rhoU.computeInteriorFaces([rho1, rho2, &u, &gamma, &faceInterpolationWeights](const Face &f) {
        Scalar flux = dot(u(f), f.outwardNorm());
        const Cell &d = flux > 0. ? f.lCell() : f.rCell();
        const Cell &a = flux > 0. ? f.rCell() : f.lCell();
        Scalar b = faceInterpolationWeights[f.id()];
        Scalar g = (1. - b) * gamma(d) + b * gamma(a);
        return ((1. - g) * rho1 + g * rho2) * u(f);
    });

    rhoU.computeBoundaryFaces([rho1, rho2, &u, &gamma](const Face &f) {
        Scalar g = gamma(f);
        return ((1. - g) * rho1 + g * rho2) * u(f);
    });
}

FiniteVolumeEquation<Scalar> cicsam::div(const VectorFiniteVolumeField &u,
                                         ScalarFiniteVolumeField &gamma,
                                         const std::vector<Scalar> &faceInterpolationWeights,
                                         Scalar theta,
                                         const CellGroup &cells)
{
    FiniteVolumeEquation<Scalar> eqn(gamma);
    const ScalarFiniteVolumeField &gamma0 = gamma.oldField(0);

    for (const Cell &cell: cells)
    {
        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar flux = dot(u(nb.face()), nb.outwardNorm());

            const Cell &donor = flux > 0. ? cell : nb.cell();
            const Cell &acceptor = flux > 0. ? nb.cell() : cell;

            //- Note, this weight is only an approximation of the correct implicit weight
            Scalar b = faceInterpolationWeights[nb.face().id()];

            eqn.add(cell, donor, theta * (1. - b) * flux);
            eqn.add(cell, acceptor, theta * b * flux);

            Scalar gammaF = (1. - b) * gamma0(donor) + b * gamma0(acceptor);
            eqn.addSource(cell, (1. - theta) * flux * gammaF);
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Scalar flux = dot(u(bd.face()), bd.outwardNorm());
            switch (gamma.boundaryType(bd.face()))
            {
            case ScalarFiniteVolumeField::FIXED:
                eqn.addSource(cell, theta * flux * gamma(bd.face()));
                eqn.addSource(cell, (1. - theta) * flux * gamma0(bd.face()));
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                eqn.add(cell, cell, theta * flux);
                eqn.addSource(cell, (1. - theta) * flux * gamma0(bd.face()));
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

FiniteVolumeEquation<Scalar> cicsam::div2e(const VectorFiniteVolumeField &u,
                                           ScalarFiniteVolumeField &gamma,
                                           const std::vector<Scalar> &faceInterpolationWeights0,
                                           const std::vector<Scalar> &faceInterpolationWeights1,
                                           Scalar theta,
                                           const CellGroup &cells)
{
    FiniteVolumeEquation<Scalar> eqn(gamma);

    for (const Cell &cell: cells)
    {
        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar flux0 = dot(u(nb.face()), nb.outwardNorm());
            Scalar flux1 = dot(u.oldField(0)(nb.face()), nb.outwardNorm());

            const Cell &donor0 = flux0 > 0. ? cell : nb.cell();
            const Cell &acceptor0 = flux0 > 0. ? nb.cell() : cell;
            const Cell &donor1 = flux1 > 0. ? cell : nb.cell();
            const Cell &acceptor1 = flux1 > 0. ? nb.cell() : cell;

            //- Note, this weight is only an approximation of the correct implicit weight
            Scalar b0 = faceInterpolationWeights0[nb.face().id()];
            Scalar b1 = faceInterpolationWeights1[nb.face().id()];

            Scalar gammaF0 = (1. - b0) * gamma(donor0) + b0 * gamma(acceptor0);
            eqn.addSource(cell, theta * flux0 * gammaF0);

            Scalar gammaF1 = (1. - b1) * gamma.oldField(0)(donor1) + b1 * gamma.oldField(0)(acceptor1);

            eqn.addSource(cell, (1 - theta) * flux1 * gammaF1);
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
                                         ScalarFiniteVolumeField &gamma,
                                         const std::vector<Scalar> &faceInterpolationWeights,
                                         Scalar theta)
{
    return div(u, gamma, faceInterpolationWeights, theta, gamma.cells());
}

FiniteVolumeEquation<Scalar> cicsam::div2e(const VectorFiniteVolumeField &u,
                                           ScalarFiniteVolumeField &gamma,
                                           const std::vector<Scalar> &faceInterpolationWeights0,
                                           const std::vector<Scalar> &faceInterpolationWeights1,
                                           Scalar theta)
{
    return div2e(u, gamma, faceInterpolationWeights0, faceInterpolationWeights1, theta, gamma.cells());
}

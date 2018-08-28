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

Scalar cicsam::beta(const VectorFiniteVolumeField &u,
                    const ScalarFiniteVolumeField &gamma,
                    const VectorFiniteVolumeField &gradGamma,
                    Scalar timeStep,
                    Scalar k,
                    const Face &f)
{
    Vector2D sf = f.outwardNorm(f.lCell().centroid());
    Scalar flux = dot(u(f), sf);
    const Cell &donor = flux >= 0. ? f.lCell() : f.rCell();
    const Cell &acceptor = flux >= 0. ? f.rCell() : f.lCell();
    const Cell &upwind = gamma.grid()->globalCells().nearestItem(2. * donor.centroid() - acceptor.centroid());

    Vector2D rc = acceptor.centroid() - donor.centroid();

    Scalar gammaD = clamp(gamma(donor), 0., 1.);
    Scalar gammaA = clamp(gamma(acceptor), 0., 1.);
    //Scalar gammaU = clamp(gammaA - 2. * dot(rc, gradGamma(donor)), 0., 1.);
    Scalar gammaU = clamp(gamma(upwind), 0., 1.);

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
    return std::isnan(betaFace) ? 0. : clamp(betaFace, 0., 1.);
}

void cicsam::computeMomentumFlux(Scalar rho1,
                                 Scalar rho2,
                                 const VectorFiniteVolumeField &u,
                                 const ScalarFiniteVolumeField &gamma,
                                 const VectorFiniteVolumeField &gradGamma,
                                 Scalar timeStep,
                                 VectorFiniteVolumeField &rhoU)
{
    const VectorFiniteVolumeField &u0 = u.oldField(0);
    const VectorFiniteVolumeField &u1 = u.oldField(1);
    const ScalarFiniteVolumeField &gamma0 = gamma.oldField(0);
    const ScalarFiniteVolumeField &gamma1 = gamma.oldField(1);
    const VectorFiniteVolumeField &gradGamma0 = gradGamma.oldField(0);
    const VectorFiniteVolumeField &gradGamma1 = gradGamma.oldField(1);

    rhoU.savePreviousTimeStep(timeStep, 2);

    rhoU.oldField(0).computeInteriorFaces([rho1, rho2, &u0, &gamma0, &gradGamma0, timeStep](const Face &f) {
        Scalar flux = dot(u0(f), f.outwardNorm());
        const Cell &d = flux > 0. ? f.lCell() : f.rCell();
        const Cell &a = flux > 0. ? f.rCell() : f.lCell();
        Scalar b = cicsam::beta(u0, gamma0, gradGamma0, timeStep, 0.5, f);
        Scalar g = (1. - b) * gamma0(d) + b * gamma0(a);
        return ((1. - g) * rho1 + g * rho2) * u0(f);
    });

    rhoU.oldField(0).computeBoundaryFaces([rho1, rho2, &u0, &gamma0](const Face &f) {
        Scalar g = gamma0(f);
        return ((1. - g) * rho1 + g * rho2) * u0(f);
    });

    rhoU.oldField(1).computeInteriorFaces([rho1, rho2, &u1, &gamma1, &gradGamma1, timeStep](const Face &f) {
        Scalar flux = dot(u1(f), f.outwardNorm());
        const Cell &d = flux > 0. ? f.lCell() : f.rCell();
        const Cell &a = flux > 0. ? f.rCell() : f.lCell();
        Scalar b = cicsam::beta(u1, gamma1, gradGamma1, timeStep, 0.5, f);
        Scalar g = (1. - b) * gamma1(d) + b * gamma1(a);
        return ((1. - g) * rho1 + g * rho2) * u1(f);
    });

    rhoU.oldField(1).computeBoundaryFaces([rho1, rho2, &u1, &gamma1](const Face &f) {
        Scalar g = gamma1(f);
        return ((1. - g) * rho1 + g * rho2) * u1(f);
    });

    rhoU.computeInteriorFaces([rho1, rho2, &u, &gamma, &gradGamma, timeStep](const Face &f) {
        Scalar flux = dot(u(f), f.outwardNorm());
        const Cell &d = flux > 0. ? f.lCell() : f.rCell();
        const Cell &a = flux > 0. ? f.rCell() : f.lCell();
        Scalar b = cicsam::beta(u, gamma, gradGamma, timeStep, 0.5, f);
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
                                         const VectorFiniteVolumeField &gradGamma,
                                         Scalar timeStep,
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

            Scalar b = beta(u, gamma, gradGamma, timeStep, 0.5, nb.face());

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
                                         ScalarFiniteVolumeField &gamma,
                                         const VectorFiniteVolumeField &gradGamma,
                                         Scalar timeStep,
                                         Scalar theta)
{
    return div(u, gamma, gradGamma, timeStep, theta, gamma.cells());
}

FiniteVolumeEquation<Scalar> cicsam::div2e(const VectorFiniteVolumeField &u,
                                           ScalarFiniteVolumeField &gamma,
                                           const VectorFiniteVolumeField &gradGamma,
                                           Scalar timeStep,
                                           Scalar theta,
                                           const CellGroup &cells)
{
    FiniteVolumeEquation<Scalar> eqn(gamma);
    const VectorFiniteVolumeField &u0 = u.oldField(0);
    const VectorFiniteVolumeField &u1 = u.oldField(1);
    const ScalarFiniteVolumeField &gamma0 = gamma;
    const ScalarFiniteVolumeField &gamma1 = gamma.oldField(0);
    const VectorFiniteVolumeField &gradGamma0 = gradGamma;
    const VectorFiniteVolumeField &gradGamma1 = gradGamma.oldField(0);

    for (const Cell &cell: cells)
    {
        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar flux0 = dot(u0(nb.face()), nb.outwardNorm());
            Scalar flux1 = dot(u1(nb.face()), nb.outwardNorm());

            const Cell &donor0 = flux0 > 0. ? cell : nb.cell();
            const Cell &acceptor0 = flux0 > 0. ? nb.cell() : cell;
            const Cell &donor1 = flux1 > 0. ? cell : nb.cell();
            const Cell &acceptor1 = flux1 > 0. ? nb.cell() : cell;

            Scalar b0 = cicsam::beta(u0, gamma0, gradGamma0, timeStep, 0.5, nb.face());
            Scalar b1 = cicsam::beta(u1, gamma1, gradGamma1, timeStep, 0.5, nb.face());

            Scalar gammaF0 = (1. - b0) * gamma0(donor0) + b0 * gamma0(acceptor0);
            Scalar gammaF1 = (1. - b1) * gamma1(donor1) + b1 * gamma1(acceptor1);

            eqn.addSource(cell, theta * flux0 * gammaF0);
            eqn.addSource(cell, (1. - theta) * flux1 * gammaF1);
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Scalar flux0 = dot(u0(bd.face()), bd.outwardNorm());
            Scalar flux1 = dot(u1(bd.face()), bd.outwardNorm());

            switch (gamma.boundaryType(bd.face()))
            {
            case ScalarFiniteVolumeField::FIXED: case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                eqn.addSource(cell, theta * flux0 * gamma0(bd.face()));
                eqn.addSource(cell, (1. - theta) * flux1 * gamma1(bd.face()));
                break;

            case ScalarFiniteVolumeField::SYMMETRY:
                break;

            default:
                throw Exception("cicsam", "div2e", "unrecognized or unspecified boundary type.");
            }
        }
    }

    return eqn;
}

FiniteVolumeEquation<Scalar> cicsam::div2e(const VectorFiniteVolumeField &u,
                                           ScalarFiniteVolumeField &gamma,
                                           const VectorFiniteVolumeField &gradGamma,
                                           Scalar timeStep,
                                           Scalar theta)
{
    return div2e(u, gamma, gradGamma, timeStep, theta, gamma.cells());
}

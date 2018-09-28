#include "Math/Algorithm.h"

#include "Cicsam.h"
#include "AxisymmetricCicsam.h"

std::vector<Scalar> axi::cicsam::faceInterpolationWeights(const VectorFiniteVolumeField &u, const ScalarFiniteVolumeField &gamma, const VectorFiniteVolumeField &gradGamma, Scalar timeStep)
{
    const Scalar k = 1;

    std::vector<Scalar> beta(gamma.grid()->faces().size(), 0.);

    for (const Face &face: gamma.grid()->interiorFaces())
    {
        Vector2D sf = face.polarOutwardNorm(face.lCell().centroid());
        Scalar flux = dot(u(face), sf);

        const Cell &d = flux > 0. ? face.lCell() : face.rCell();
        const Cell &a = flux <= 0. ? face.lCell() : face.rCell();
        //const Cell &upwind = gamma.grid()->globalCells().nearestItem(2. * donor.centroid() - acceptor.centroid());

        Vector2D rc = a.centroid() - d.centroid();

        Scalar gammaD = clamp(gamma(d), 0., 1.);
        Scalar gammaA = clamp(gamma(a), 0., 1.);
        Scalar gammaU = clamp(gammaA - 2. * dot(rc, gradGamma(d)), 0., 1.);
        // Scalar gammaU = clamp(gamma(upwind), 0., 1.);

        Scalar gammaDTilde = (gammaD - gammaU) / (gammaA - gammaU);

        Scalar coD = 0.; //- Cell courant number
        for (const InteriorLink &nb: d.neighbours())
            coD += std::max(dot(u(nb.face()), nb.polarOutwardNorm()) / d.polarVolume() * timeStep, 0.);

        for (const BoundaryLink &bd: d.boundaries())
            coD += std::max(dot(u(bd.face()), bd.polarOutwardNorm()) / d.polarVolume() * timeStep, 0.);

        Scalar thetaF = std::acos(std::abs(dot(gradGamma(d).unitVec(), rc.unitVec())));
        Scalar psiF = std::min(k * (std::cos(2 * thetaF) + 1.) / 2., 1.);
        Scalar gammaFTilde = psiF * ::cicsam::hc(gammaDTilde, coD) + (1. - psiF) * ::cicsam::uq(gammaDTilde, coD);
        Scalar betaFace = (gammaFTilde - gammaDTilde) / (1. - gammaDTilde);

        //- If stencil cannot be computed, default to upwind
        beta[face.id()] = std::isnan(betaFace) ? 0. : clamp(betaFace, 0., 1.);
    }

    return beta;
}

FiniteVolumeEquation<Scalar> axi::cicsam::div(const VectorFiniteVolumeField &u, ScalarFiniteVolumeField &gamma, const std::vector<Scalar> &faceInterpolationWeights, Scalar theta, const CellGroup &cells)
{
    FiniteVolumeEquation<Scalar> eqn(gamma, 5);
    const ScalarFiniteVolumeField &gamma0 = gamma.oldField(0);

    for(const Cell &cell: gamma.cells())
    {
        for(const InteriorLink &nb: cell.neighbours())
        {
            Scalar b = faceInterpolationWeights[nb.face().id()];

            Scalar flux = dot(u(nb.face()), nb.polarOutwardNorm());

            const Cell &d = flux > 0. ? cell : nb.cell();
            const Cell &a = flux <= 0. ? cell : nb.cell();

            eqn.add(cell, d, (1. - b) * flux * theta);
            eqn.add(cell, a, b * flux * theta);
            eqn.addSource(cell, flux * ((1. - b) * gamma0(d) + b * gamma0(a)) * (1. - theta));
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            Scalar flux = dot(u(bd.face()), bd.polarOutwardNorm());

            switch(gamma.boundaryType(bd.face()))
            {
            case ScalarFiniteVolumeField::FIXED:
                eqn.addSource(cell, flux * gamma(bd.face()) * theta);
                eqn.addSource(cell, flux * gamma0(cell) * (1. - theta));
                break;
            case ScalarFiniteVolumeField::SYMMETRY:
            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                break;
            default:
                throw Exception("axi::cicsam", "div", "unrecognized boundary condition type.");
            }
        }
    }

    return eqn;
}

void axi::cicsam::computeMomentumFlux(Scalar rho1,
                                      Scalar rho2,
                                      const VectorFiniteVolumeField &u,
                                      const ScalarFiniteVolumeField &gamma,
                                      const std::vector<Scalar> &faceInterpolationWeights,
                                      VectorFiniteVolumeField &rhoU)
{
    for(const Face &f: rhoU.grid()->interiorFaces())
    {
        Scalar flux = dot(u(f), f.polarOutwardNorm(f.lCell().centroid()));

        const Cell &d = flux > 0. ? f.lCell() : f.rCell();
        const Cell &a = flux <= 0. ? f.lCell() : f.rCell();

        Scalar b = faceInterpolationWeights[f.id()];
        Scalar gf = (1. - b) * gamma(d) + b * gamma(a);

        rhoU(f) = ((1. - gf) * rho1 + gf * rho2) * u(f);
    }

    for(const Face &f: rhoU.grid()->boundaryFaces())
    {
        Scalar gf = gamma(f);
        rhoU(f) = ((1. - gf) * rho1 + gf * rho2) * u(f);
    }
}

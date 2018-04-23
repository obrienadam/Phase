#include "Hric.h"
#include "Math/Algorithm.h"
#include "Cicsam.h"

ScalarFiniteVolumeField hric::beta(const VectorFiniteVolumeField &u,
                                   const VectorFiniteVolumeField &gradGamma,
                                   const ScalarFiniteVolumeField &gamma,
                                   Scalar timeStep)
{
    ScalarFiniteVolumeField beta(gamma.grid(), "beta");

    for(const Face& face: gamma.grid()->interiorFaces())
    {
        Vector2D sf = face.outwardNorm(face.lCell().centroid());
        Scalar flux = dot(u(face), sf);
        const Cell& donor = flux >= 0. ? face.lCell() : face.rCell();
        const Cell& acceptor = flux >= 0. ? face.rCell() : face.lCell();
        Vector2D rc = acceptor.centroid() - donor.centroid();

        Scalar gammaD = clamp(gamma(donor), 0., 1.);
        Scalar gammaA = clamp(gamma(acceptor), 0., 1.);
        Scalar gammaU = clamp(gammaA - 2.*dot(rc, gradGamma(donor)), 0., 1.);
        Scalar gammaDTilde = (gammaD - gammaU) / (gammaA - gammaU);

        Scalar coD = 0.; //- Cell courant number
        for(const InteriorLink& nb: donor.neighbours())
            coD += std::max(dot(u(nb.face()), nb.outwardNorm()) / donor.volume() * timeStep, 0.);

        for(const BoundaryLink& bd: donor.boundaries())
            coD += std::max(dot(u(bd.face()), bd.outwardNorm()) / donor.volume() * timeStep, 0.);

        Scalar gammaFTilde = gammaDTilde < 0. || gammaDTilde > 1. ? gammaDTilde:
                             0. <= gammaDTilde && gammaDTilde < 0.5 ? 2. * gammaDTilde: 1.;

        Scalar lambdaF = std::sqrt(std::abs(dot(gradGamma(donor).unitVec(), rc.unitVec())));

        gammaFTilde = lambdaF * gammaFTilde + (1. - lambdaF) * gammaDTilde;
        gammaFTilde = coD < 0.3 ? gammaFTilde:
                      coD > 0.7 ? gammaDTilde:
                      gammaDTilde + (gammaFTilde - gammaDTilde) * (0.7 - coD) / (0.7 - 0.3);

        Scalar betaFace = (gammaFTilde - gammaDTilde) / (1. - gammaDTilde);

        //- If stencil cannot be computed, default to upwind
        beta(face) = std::isnan(betaFace) ? 0.: clamp(betaFace, 0., 1.);
    }

    return beta;
}

FiniteVolumeEquation<Scalar> hric::div(const VectorFiniteVolumeField &u,
                           const ScalarFiniteVolumeField& beta,
                           ScalarFiniteVolumeField &gamma,
                           Scalar theta)
{
    return cicsam::div(u, beta, gamma, theta);
}
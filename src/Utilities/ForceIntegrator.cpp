#include <boost/algorithm/string.hpp>

#include "ForceIntegrator.h"

//- Public static functions

std::vector<ForceIntegrator> ForceIntegrator::initForceIntegrators(const Input &input,
                                                                   const ScalarFiniteVolumeField &p,
                                                                   const ScalarFiniteVolumeField &rho,
                                                                   const ScalarFiniteVolumeField &mu,
                                                                   const VectorFiniteVolumeField &u)
{
    std::vector<ForceIntegrator> forceIntegrators;

    boost::optional<std::string> forceIntegratorInput = input.caseInput().get_optional<std::string>("Integrators.ForceIntegrators");

    if(!forceIntegratorInput)
        return forceIntegrators;

    std::string patchInput = input.caseInput().get<std::string>("Integrators.ForceIntegrators.patches");
    std::vector<std::string> patchNames;

    boost::split(patchNames, patchInput, boost::is_any_of(", "), boost::token_compress_on);

    for(const std::string& patchName: patchNames)
    {
        printf("Initializing a force integrator on patch \"%s\".\n", patchName.c_str());

        forceIntegrators.push_back(
                    ForceIntegrator(p.grid().patch(patchName), p, rho, mu, u)
                    );
    }

    return forceIntegrators;
}

// Constructors
ForceIntegrator::ForceIntegrator(const Patch &patch,
                                 const ScalarFiniteVolumeField &p,
                                 const ScalarFiniteVolumeField &rho,
                                 const ScalarFiniteVolumeField &mu,
                                 const VectorFiniteVolumeField &u)
    :
      patch_(patch),
      p_(p),
      rho_(rho),
      mu_(mu),
      u_(u)
{

}

Vector2D ForceIntegrator::integrate()
{
    Vector2D fn = Vector2D(0., 0.), ft = Vector2D(0., 0.);

    const auto sqr = [](Scalar x) { return x*x; };

    for(const Face &face: patch_)
    {
        const Cell &cell = face.lCell();
        const Vector2D sfn = face.outwardNorm(cell.centroid());
        const Vector2D sft = -sfn.tangentVec();

        fn += (p_(face) + 0.5*rho_(face)*sqr(dot(u_(face), sfn))/sfn.magSqr())*sfn;
        ft += mu_(cell)*dot(u_(cell) - u_(face), sft.unitVec())*sft;
    }

    printf("Net normal force on patch \"%s\": %s\n", patch_.name().c_str(), to_string(fn).c_str());
    printf("Net tangential force on patch \"%s\": %s\n", patch_.name().c_str(), to_string(ft).c_str());
    printf("Net force on patch \"%s\": %s\n", patch_.name().c_str(), to_string(fn + ft).c_str());

    data_.push_back(fn + ft);

    return data_.back();
}

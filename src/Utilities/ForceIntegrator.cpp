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
    data_.push_back(computeForce(patch_, p_, rho_, mu_, u_));
    return data_.back();
}

Vector2D computeForce(const FaceGroup& patch,
                      Scalar rho,
                      Scalar mu,
                      const ScalarFiniteVolumeField& p,
                      const VectorFiniteVolumeField &u)
{
    Vector2D fn = Vector2D(0., 0.), ft = Vector2D(0., 0.);

    for(const Face &face: patch)
    {
        const Cell &cell = face.lCell();
        Vector2D sf = face.outwardNorm(cell.centroid());
        Vector2D tf = sf.tangentVec();
        Vector2D sfn = sf.unitVec();
        Scalar l = (cell.centroid() - face.centroid()).mag();

        fn -= (p(face) + 0.5*rho*pow(dot(u(face), sfn), 2))*sf;
        ft += mu*dot(u(cell) - u(face), tf)*tf/(tf.magSqr()*l);
    }

    printf("Net normal force on patch \"%s\": %s\n", patch.name().c_str(), std::to_string(fn).c_str());
    printf("Net tangential force on patch \"%s\": %s\n", patch.name().c_str(), std::to_string(ft).c_str());
    printf("Net force on patch \"%s\": %s\n", patch.name().c_str(), std::to_string(fn + ft).c_str());

    return fn + ft;
}

Vector2D computeForce(const Patch& patch,
                      const ScalarFiniteVolumeField& p,
                      const ScalarFiniteVolumeField& rho,
                      const ScalarFiniteVolumeField& mu,
                      const VectorFiniteVolumeField &u)
{
    Vector2D fn = Vector2D(0., 0.), ft = Vector2D(0., 0.);

    for(const Face &face: patch)
    {
        const Cell &cell = face.lCell();
        Vector2D sf = face.outwardNorm(cell.centroid());
        Vector2D tf = sf.tangentVec();
        Vector2D sfn = sf.unitVec();
        Scalar l = (cell.centroid() - face.centroid()).mag();

        fn -= (p(face) + 0.5*rho(face)*pow(dot(u(face), sfn), 2))*sf;
        ft += 0.5*(mu(cell) + mu(face))*dot(u(cell) - u(face), tf)*tf/(tf.magSqr()*l);
    }

    printf("Net normal force on patch \"%s\": %s\n", patch.name().c_str(), std::to_string(fn).c_str());
    printf("Net tangential force on patch \"%s\": %s\n", patch.name().c_str(), std::to_string(ft).c_str());
    printf("Net force on patch \"%s\": %s\n", patch.name().c_str(), std::to_string(fn + ft).c_str());

    return fn + ft;
}

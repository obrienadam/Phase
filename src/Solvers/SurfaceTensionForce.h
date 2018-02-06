#ifndef SURFACE_TENSION_FORCE_H
#define SURFACE_TENSION_FORCE_H

#include "Solver.h"
#include "ImmersedBoundary.h"
#include "VectorFiniteVolumeField.h"
#include "ScalarGradient.h"

class SurfaceTensionForce : public VectorFiniteVolumeField
{
public:

    //- Constructor
    SurfaceTensionForce(const Input &input,
                        const std::weak_ptr<ImmersedBoundary> &ib,
                        ScalarFiniteVolumeField &gamma,
                        const ScalarFiniteVolumeField &rho,
                        const ScalarFiniteVolumeField &mu,
                        const VectorFiniteVolumeField &u,
                        const ScalarGradient &gradGamma);

    //- Internal field pointers
    const std::shared_ptr<ScalarFiniteVolumeField> &kappa() const
    { return kappa_; }

    const std::shared_ptr<ScalarFiniteVolumeField> &gammaTilde() const
    { return gammaTilde_; }

    const std::shared_ptr<VectorFiniteVolumeField> &n() const
    { return n_; }

    const std::shared_ptr<ScalarGradient> &gradGammaTilde() const
    { return gradGammaTilde_; }

    //- Return properties
    Scalar sigma() const
    { return sigma_; }

    Scalar theta(const ImmersedBoundaryObject &ibObj) const
    {
        auto it = ibContactAngles_.find(ibObj.id());
        return it != ibContactAngles_.end() ? it->second : M_PI_2;
    }

    Scalar theta(const Patch &patch) const
    {
        auto it = patchContactAngles_.find(patch.id());
        return it != patchContactAngles_.end() ? it->second : M_PI_2;
    }

    //- Compute
    virtual void computeFaces() = 0;

    virtual void compute() = 0;

    virtual void computeInterfaceNormals();

    Vector2D contactLineNormal(const Cell &cell, const Point2D &pt,
                               const ImmersedBoundaryObject &ibObj) const;

    Vector2D contactLineNormal(const Cell &cell,
                               const ImmersedBoundaryObject &ibObj) const
    { return contactLineNormal(cell, ibObj.nearestIntersect(cell.centroid()), ibObj); }

    void smoothGammaField();

    //- Misc special gamma boundary equations
    virtual Equation<Scalar> contactLineBcs();

protected:

    Scalar sigma_, kernelWidth_, eps_ = 1e-8;

    std::unordered_map<Label, Scalar> ibContactAngles_, patchContactAngles_;

    const ScalarFiniteVolumeField &rho_;
    const ScalarFiniteVolumeField &mu_;
    const VectorFiniteVolumeField &u_;
    const ScalarGradient &gradGamma_;
    ScalarFiniteVolumeField &gamma_;

    std::weak_ptr<ImmersedBoundary> ib_;

    //- Fields, can share ownership
    std::shared_ptr<ScalarFiniteVolumeField> kappa_, gammaTilde_;
    std::shared_ptr<VectorFiniteVolumeField> n_;
    std::shared_ptr<ScalarGradient> gradGammaTilde_;
};

#endif

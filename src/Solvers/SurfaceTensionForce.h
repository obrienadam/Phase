#ifndef SURFACE_TENSION_FORCE_H
#define SURFACE_TENSION_FORCE_H

#include "Solver.h"
#include "ImmersedBoundary.h"
#include "VectorFiniteVolumeField.h"
#include "ScalarGradient.h"

class SurfaceTensionForce: public VectorFiniteVolumeField
{
public:

    //- Constructor
    SurfaceTensionForce(const Input &input,
                        const ImmersedBoundary &ib,
                        const ScalarFiniteVolumeField& gamma,
                        const ScalarFiniteVolumeField& rho,
                        const ScalarFiniteVolumeField& mu,
                        const VectorFiniteVolumeField& u,
                        const ScalarGradient& gradGamma);

    //- Check if a particular patch should have contact lines enforced
    bool isContactLinePatch(const Patch &patch) const
    { return patchContactAngles_.find(patch.id()) != patchContactAngles_.end(); }

    bool isContactLineIbObj(const ImmersedBoundaryObject& ibObj) const
    { return ibContactAngles_.find(ibObj.id()) != ibContactAngles_.end(); }

    //- References to interally stored fields
    virtual void registerSubFields(Solver& solver) = 0;

    const ScalarFiniteVolumeField &kappa() const
    { return *kappa_; }

    const VectorFiniteVolumeField &n() const
    { return *n_; }

    //- Return properties
    Scalar sigma() const
    { return sigma_; }

    Scalar theta(const ImmersedBoundaryObject& ibObj) const
    {
        auto it = ibContactAngles_.find(ibObj.id());
        return it != ibContactAngles_.end() ? it->second: M_PI_2;
    }

    Scalar theta(const Patch& patch) const
    {
        auto it = patchContactAngles_.find(patch.id());
        return it != patchContactAngles_.end() ? it->second: M_PI_2;
    }

    //- Compute
    virtual void compute() = 0;

    virtual void computeInterfaceNormals();

    Vector2D contactLineNormal(const Face& face);

    Vector2D contactLineNormal(const Cell& lCell,
                               const Cell& rCell,
                               const ImmersedBoundaryObject& ibObj);

protected:

    Scalar sigma_, kernelWidth_;

    std::unordered_map<Label, Scalar> ibContactAngles_;
    std::unordered_map<Label, Scalar> patchContactAngles_;

    const ScalarFiniteVolumeField &gamma_;
    const ScalarFiniteVolumeField &rho_;
    const ScalarFiniteVolumeField &mu_;
    const VectorFiniteVolumeField &u_;
    const ScalarGradient &gradGamma_;

    const ImmersedBoundary &ib_;

    std::shared_ptr<ScalarFiniteVolumeField> kappa_, gammaTilde_;
    std::shared_ptr<ScalarGradient> gradGammaTilde_;
    std::shared_ptr<VectorFiniteVolumeField> n_;
};

#endif

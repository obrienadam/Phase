#ifndef SURFACE_TENSION_FORCE_H
#define SURFACE_TENSION_FORCE_H

#include "VectorFiniteVolumeField.h"
#include "ImmersedBoundary.h"
#include "Solver.h"

class SurfaceTensionForce: public VectorFiniteVolumeField
{
public:

    //- Constructor
    SurfaceTensionForce(const Input &input,
                        const ScalarFiniteVolumeField& gamma,
                        const ScalarFiniteVolumeField& rho,
                        const ScalarFiniteVolumeField& mu,
                        const VectorFiniteVolumeField& u,
                        VectorFiniteVolumeField& gradGamma);

    //- Compute
    virtual void compute() = 0;

    Scalar theta(const Cell& cell);

    Vector2D computeContactLineNormal(const Vector2D &gradGamma, const Vector2D &wallNormal, const Vector2D &vel,
                                      Scalar theta) const;

    Vector2D computeContactLineNormal(const Vector2D &gradGamma, const Vector2D &wallNormal, const Vector2D &vel) const;

    void computeGhostCellVals(const ImmersedBoundaryObject& ibObj, Scalar theta);

    //- Check if a particular patch should have contact lines enforced
    bool isContactLinePatch(const Patch &patch) const;

    //- References to interally stored fields
    virtual const ScalarFiniteVolumeField &gammaTilde() const = 0;

    virtual const VectorFiniteVolumeField &gradGammaTilde() const = 0;

    virtual void registerSubFields(Solver& solver) = 0;

    const ScalarFiniteVolumeField &kappa() const
    { return *kappa_; }

    const VectorFiniteVolumeField &n() const
    { return *n_; }

    //- Return properties
    Scalar theta() const
    { return thetaAdv_; }

    Scalar thetaAdv() const
    { return thetaAdv_; }

    Scalar thetaRec() const
    { return thetaRec_; }

    Scalar ibTheta(const ImmersedBoundaryObject &ibObj) const;

    Scalar sigma() const
    { return sigma_; }

protected:

    Scalar sigma_, thetaAdv_, thetaRec_;
    std::map<Label, Scalar> ibContactAngles_;

    const ScalarFiniteVolumeField &gamma_;
    const ScalarFiniteVolumeField &rho_;
    const ScalarFiniteVolumeField &mu_;
    const VectorFiniteVolumeField &u_;
    VectorFiniteVolumeField &gradGamma_;

    std::shared_ptr<ScalarFiniteVolumeField> kappa_;
    std::shared_ptr<VectorFiniteVolumeField> n_;

    std::vector<Ref<const Patch> > contactAnglePatches_;
};

#endif

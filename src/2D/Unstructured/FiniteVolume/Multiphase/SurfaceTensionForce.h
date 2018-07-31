#ifndef PHASE_SURFACE_TENSION_FORCE_H
#define PHASE_SURFACE_TENSION_FORCE_H

#include "FiniteVolume/Field/VectorFiniteVolumeField.h"
#include "FiniteVolume/Field/ScalarGradient.h"

class SurfaceTensionForce
{
public:

    //- Constructor
    SurfaceTensionForce(const Input &input,
                        const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                        const std::shared_ptr<const CellGroup> &fluid);

    virtual void setCellGroup(const std::shared_ptr<const CellGroup> &fluid);

    //- Internal field pointers
    const std::shared_ptr<VectorFiniteVolumeField> &fst() const
    { return fst_; }

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

    Scalar theta(const FaceGroup &patch) const;

    //- Compute
    virtual void computeFaceInterfaceForces(const ScalarFiniteVolumeField &gamma, const ScalarGradient &gradGamma) = 0;

    virtual void computeInterfaceForces(const ScalarFiniteVolumeField &gamma, const ScalarGradient &gradGamma) = 0;

    virtual void computeInterfaceNormals();

    void smoothGammaField(const ScalarFiniteVolumeField &gamma);

protected:

    std::shared_ptr<const FiniteVolumeGrid2D> grid_;

    std::shared_ptr<const CellGroup> fluid_;

    Scalar sigma_, kernelWidth_, eps_ = 1e-8;

    std::unordered_map<std::string, Scalar> patchContactAngles_;

    //- Fields, can share ownership
    std::shared_ptr<VectorFiniteVolumeField> fst_;

    std::shared_ptr<ScalarFiniteVolumeField> kappa_, gammaTilde_;

    std::shared_ptr<VectorFiniteVolumeField> n_;

    std::shared_ptr<ScalarGradient> gradGammaTilde_;
};

#endif

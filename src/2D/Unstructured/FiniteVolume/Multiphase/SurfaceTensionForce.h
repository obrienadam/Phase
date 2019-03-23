#ifndef PHASE_SURFACE_TENSION_FORCE_H
#define PHASE_SURFACE_TENSION_FORCE_H

#include "FiniteVolume/Field/VectorFiniteVolumeField.h"
#include "FiniteVolume/Field/ScalarGradient.h"

class SurfaceTensionForce
{
public:

    class SmoothingKernel
    {
    public:

        enum Type{PESKIN, POW_6, POW_8};

        SmoothingKernel(const Cell& cell, Scalar eps, Type type = POW_8);

        void setAxisymmetric(bool axisymmetric);

        const Cell &cell() const
        { return cell_; }

        Scalar eval(const ScalarFiniteVolumeField &phi) const;

    private:

        //        Scalar kernel(Scalar r) const
        //        { return r < eps_ ? std::cos(M_PI * r / eps_) + 1. : 0.; }

        Scalar kcos(Scalar x) const
        { return x < eps_ ?  eps_ * (1. + std::cos(M_PI * x / eps_)) : 0.; }

        Scalar kernel(Vector2D dx, Type type) const;

        bool axisymmetric_ = false;

        Type type_;

        Scalar eps_, A_;

        const Cell &cell_;

        std::vector<Ref<const Cell>> kCells_;
    };

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

    //- Helpers
    Scalar dynamicContactAngle(Scalar thetaApp, Scalar Ca, Scalar delta, Scalar K = 0.02) const;

    //- Compute
    virtual void computeFaceInterfaceForces(const ScalarFiniteVolumeField &gamma, const ScalarGradient &gradGamma) = 0;

    virtual void computeInterfaceForces(const ScalarFiniteVolumeField &gamma, const ScalarGradient &gradGamma) = 0;

    virtual void computeInterfaceNormals();

    void setAxisymmetric(bool axisymmetric);

    void smoothGammaField(const ScalarFiniteVolumeField &gamma);

protected:

    static SmoothingKernel::Type getKernelType(std::string type);

    std::shared_ptr<const FiniteVolumeGrid2D> grid_;

    std::shared_ptr<const CellGroup> fluid_;

    Scalar sigma_, kernelWidth_, eps_ = 1e-8;

    SmoothingKernel::Type kernelType_;

    Scalar minTheta_, maxTheta_;

    std::unordered_map<std::string, Scalar> patchContactAngles_;

    std::vector<SmoothingKernel> kernels_;

    //- Fields, can share ownership
    std::shared_ptr<VectorFiniteVolumeField> fst_;

    std::shared_ptr<ScalarFiniteVolumeField> kappa_, gammaTilde_;

    std::shared_ptr<VectorFiniteVolumeField> n_;

    std::shared_ptr<ScalarGradient> gradGammaTilde_;
};

#endif

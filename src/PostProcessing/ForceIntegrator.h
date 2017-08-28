#ifndef FORCE_INTEGRATOR_H
#define FORCE_INTEGRATOR_H

#include "Patch.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

class ForceIntegrator
{
public:

    ForceIntegrator(const Patch& patch,
                    const VectorFiniteVolumeField& u,
                    const ScalarFiniteVolumeField& rho,
                    const ScalarFiniteVolumeField& mu,
                    const ScalarFiniteVolumeField& p);

    void compute(Scalar time);

    void write() const;

private:

    const Patch& patch_;
    const VectorFiniteVolumeField& u_;
    const ScalarFiniteVolumeField &rho_, &mu_, &p_;

    std::vector<Scalar> time_;
    std::vector<Vector2D> force_;
};

#endif

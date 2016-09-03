#ifndef FORCE_INTEGRATOR_H
#define FORCE_INTEGRATOR_H

#include "Patch.h"
#include "Vector2D.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"
#include "ImmersedBoundaryObject.h"

class ForceIntegrator
{
public:

    static std::vector<ForceIntegrator> initForceIntegrators(const Input &input,
                                                             const ScalarFiniteVolumeField& p,
                                                             const ScalarFiniteVolumeField& rho,
                                                             const ScalarFiniteVolumeField& mu,
                                                             const VectorFiniteVolumeField &u);

    ForceIntegrator(const Patch &patch,
                    const ScalarFiniteVolumeField& p,
                    const ScalarFiniteVolumeField& rho,
                    const ScalarFiniteVolumeField& mu,
                    const VectorFiniteVolumeField& u);

    Vector2D integrate();

private:

    const Patch &patch_;
    const ScalarFiniteVolumeField &p_, &rho_, &mu_;
    const VectorFiniteVolumeField &u_;

    std::vector<Vector2D> data_;
};

#endif

#ifndef VOLUME_INTEGRATOR_H
#define VOLUME_INTEGRATOR_H

#include "PostProcessingObject.h"

class VolumeIntegrator : public PostProcessingObject
{
public:

    VolumeIntegrator(const Solver &solver, const std::string &fieldName);

    void compute(Scalar time);

private:

    const ScalarFiniteVolumeField& field_;

};


#endif //PHASE_VOLUMEINTEGRATOR_H

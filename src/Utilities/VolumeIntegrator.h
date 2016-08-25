#ifndef VOLUME_INTEGRATOR_H
#define VOLUME_INTEGRATOR_H

#include "ScalarFiniteVolumeField.h"

class Solver;

class VolumeIntegrator
{
public:

    static std::vector<VolumeIntegrator> initVolumeIntegrators(const Input& input, const Solver& solver);

    VolumeIntegrator(const ScalarFiniteVolumeField& field, const CellGroup &cellGroup);

    Scalar integrate() const;

private:

    const ScalarFiniteVolumeField &field_;
    const CellGroup &cellGroup_;

};

#endif

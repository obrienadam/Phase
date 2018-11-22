#ifndef PHASE_POST_PROCESSING_H
#define PHASE_POST_PROCESSING_H

#include "System/PostProcessingInterface.h"
#include "System/Input.h"
#include "Solvers/Solver.h"
#include "FiniteVolume/ImmersedBoundary/ImmersedBoundary.h"

#include "Viewer.h"

class PostProcessing : public PostProcessingInterface
{
public:

    PostProcessing(const Input& input, const Solver& solver);

    void initIbPostProcessingObjects(const Input &input, const Solver &solver);

    void compute(Scalar time, bool force = false) override;

protected:

    int iter_, fileWriteFrequency_;

    std::unique_ptr<Viewer> viewer_;
};


#endif //PHASE_POSTPROCESSING_H

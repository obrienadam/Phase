#ifndef PHASE_POST_PROCESSING_H
#define PHASE_POST_PROCESSING_H

#include "System/PostProcessingInterface.h"
#include "System/Input.h"
#include "Solvers/Solver.h"

#include "CgnsViewer.h"

class PostProcessing : public PostProcessingInterface
{
public:

    PostProcessing(const Input& input, const Solver& solver);

    void compute(Scalar time);

protected:

    int iter_, fileWriteFrequency_;

    CgnsViewer viewer_;
};


#endif //PHASE_POSTPROCESSING_H

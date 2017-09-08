#ifndef POST_PROCESSING_H
#define POST_PROCESSING_H

#include <memory>

#include "PostProcessingObject.h"
#include "Input.h"

class PostProcessing
{

public:

    PostProcessing(const Input &input, const Solver &solver);

    void compute(Scalar time) const;

private:

    std::vector<std::shared_ptr<PostProcessingObject>> postProcessingObjs_;
};

#endif
#ifndef POST_PROCESSING_H
#define POST_PROCESSING_H

#include <memory>
#include <vector>

#include "Types.h"
#include "Input.h"

class Solver;

class PostProcessing
{
public:

    static std::vector<std::shared_ptr<PostProcessing>> initPostProcessing(const Input& input, const Solver& solver);

    PostProcessing(const Solver& solver);

    virtual void compute(Scalar time) = 0;

protected:

    const Solver& solver_;
};


#endif

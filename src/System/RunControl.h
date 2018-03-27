#ifndef RUN_CONTROL_H
#define RUN_CONTROL_H

#include "Input.h"
#include "SolverInterface.h"
#include "PostProcessingInterface.h"
#include "Timer.h"

class RunControl
{
public:

    void run(const Input &input,
             SolverInterface &solver,
             PostProcessingInterface &postProcessing);

private:
    Timer time_;
};

#endif

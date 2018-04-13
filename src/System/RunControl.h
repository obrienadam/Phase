#ifndef RUN_CONTROL_H
#define RUN_CONTROL_H

#include "CommandLine.h"
#include "Input.h"
#include "SolverInterface.h"
#include "PostProcessingInterface.h"
#include "Timer.h"

class RunControl
{
public:

    void run(const CommandLine &cl,
             const Input &input,
             SolverInterface &solver,
             PostProcessingInterface &postProcessing);

private:
    Timer time_;
};

#endif

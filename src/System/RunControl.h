#ifndef PHASE_RUN_CONTROL_H
#define PHASE_RUN_CONTROL_H

#include "CommandLine.h"
#include "Input.h"
#include "PostProcessingInterface.h"
#include "SolverInterface.h"
#include "Timer.h"

class RunControl {
public:
  void run(const CommandLine &cl, const Input &input, SolverInterface &solver,
           PostProcessingInterface &postProcessing);

private:
  Timer time_;
};

#endif

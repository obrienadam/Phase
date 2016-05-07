#ifndef RUN_CONTROL_H
#define RUN_CONTROL_H

#include "Input.h"
#include "Solver.h"
#include "Viewer.h"

class RunControl
{
public:

    void run(const Input& input,
             Solver & solver,
             Viewer &viewer);
};

#endif

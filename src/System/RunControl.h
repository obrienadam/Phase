#ifndef RUN_CONTROL_H
#define RUN_CONTROL_H

#include "Input.h"
#include "Communicator.h"
#include "Solver.h"
#include "CgnsViewer.h"
#include "Time.h"

class RunControl
{
public:

    void run(const Input& input,
             const Communicator& comm,
             Solver & solver,
             Viewer &viewer);

private:
    Time time_;
};

#endif

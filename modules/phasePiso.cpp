#include <iostream>

#include "Input.h"
#include "Communicator.h"
#include "CommandLine.h"
#include "ConstructGrid.h"
#include "Piso.h"
#include "Viewer.h"
#include "RunControl.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    Input input;
    Communicator comm;
    CommandLine(argc, argv);

    input.parseInputFile();

    shared_ptr<FiniteVolumeGrid2D> gridPtr(constructGrid(input));

    Piso solver(input, comm, *gridPtr);

    Viewer viewer(solver, input);
    RunControl runControl;

    runControl.run(input, solver, viewer);

    return 0;
}


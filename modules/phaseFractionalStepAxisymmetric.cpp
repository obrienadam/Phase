#include <iostream>

#include "Input.h"
#include "Communicator.h"
#include "CommandLine.h"
#include "FractionalStepAxisymmetric.h"
#include "CgnsViewer.h"
#include "RunControl.h"

int main(int argc, char *argv[])
{
    using namespace std;

    Communicator::init(argc, argv);

    Input input;
    CommandLine(argc, argv);

    input.parseInputFile();

    FractionalStepAxisymmetric solver(input);
    CgnsViewer viewer(input, solver);

    RunControl runControl;
    runControl.run(input, solver, viewer);

    Communicator::finalize();

    return 0;
}


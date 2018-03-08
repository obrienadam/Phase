#include <iostream>

#include "Input.h"
#include "CommandLine.h"
#include "SolverFactory.h"
#include "CgnsViewer.h"
#include "RunControl.h"

int main(int argc, char *argv[])
{
    using namespace std;

    Communicator::init(argc, argv);

    Input input;
    CommandLine(argc, argv);

    input.parseInputFile();

    shared_ptr<Solver> solver = SolverFactory::create(SolverFactory::FRACTIONAL_STEP, input);

    CgnsViewer viewer(input, *solver);
    RunControl runControl;

    runControl.run(input, *solver, viewer);

    Communicator::finalize();

    return 0;
}


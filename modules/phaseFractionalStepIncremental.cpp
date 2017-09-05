#include <iostream>

#include "Input.h"
#include "Communicator.h"
#include "CommandLine.h"
#include "ConstructGrid.h"
#include "FractionalStepIncremental.h"
#include "CgnsViewer.h"
#include "RunControl.h"

int main(int argc, char *argv[])
{
    using namespace std;

    Communicator::init(argc, argv);

    Input input;
    CommandLine(argc, argv);

    input.parseInputFile();

    shared_ptr<FiniteVolumeGrid2D> grid = constructGrid(input, std::make_shared<Communicator>());

    FractionalStepIncremental solver(input, grid);
    CgnsViewer viewer(input, solver);

    RunControl runControl;
    runControl.run(input, solver, viewer);

    Communicator::finalize();

    return 0;
}


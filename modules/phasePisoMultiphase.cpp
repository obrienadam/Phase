#include <iostream>

#include "Input.h"
#include "Communicator.h"
#include "CommandLine.h"
#include "ConstructGrid.h"
#include "PisoMultiphase.h"
#include "CgnsViewer.h"
#include "RunControl.h"

int main(int argc, char *argv[])
{
    using namespace std;

    Communicator::init(argc, argv);

    Input input;
    CommandLine(argc, argv);

    input.parseInputFile();

    shared_ptr<FiniteVolumeGrid2D> grid(constructGrid(input));
    grid->partition(input, std::make_shared<Communicator>());

    PisoMultiphase solver(input, grid);
    CgnsViewer viewer(input, solver);

    RunControl runControl;

    runControl.run(input, solver, viewer);

    Communicator::finalize();

    return 0;
}


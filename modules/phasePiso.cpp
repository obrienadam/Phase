#include <iostream>

#include "Input.h"
#include "Communicator.h"
#include "CommandLine.h"
#include "ConstructGrid.h"
#include "Piso.h"
#include "CgnsViewer.h"
#include "RunControl.h"

int main(int argc, char *argv[])
{
    using namespace std;

    Communicator::init(argc, argv);

    Input input;
    Communicator comm;
    CommandLine(argc, argv);

    input.parseInputFile();

    shared_ptr<FiniteVolumeGrid2D> grid(constructGrid(input));
    grid->partition(input, comm);

    Piso solver(input, comm, grid);
    CgnsViewer viewer(input, comm, solver);

    RunControl runControl;
    runControl.run(input, comm, solver, viewer);

    Communicator::finalize();

    return 0;
}


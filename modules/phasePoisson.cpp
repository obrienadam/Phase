#include <iostream>

#include "Input.h"
#include "CommandLine.h"
#include "ConstructGrid.h"
#include "Viewer.h"
#include "Poisson.h"

int main(int argc, char* argv[])
{
    using namespace std;

    Communicator::init(argc, argv);

    Input input;
    Communicator comm;
    CommandLine(argc, argv);

    input.parseInputFile();

    shared_ptr<FiniteVolumeGrid2D> gridPtr(constructGrid(input));
    gridPtr->partition(comm);


    Poisson solver(input, comm, *gridPtr);
    Viewer viewer(input, comm, solver);

    //viewer.write(1., comm);

    Communicator::finalize();

    return 0;
}

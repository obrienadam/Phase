#include <iostream>

#include "Input.h"
#include "CommandLine.h"
#include "ConstructGrid.h"
#include "CgnsViewer.h"
#include "Poisson.h"
#include "Time.h"

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
    CgnsViewer viewer(input, comm, solver);

    Time time;

    time.start();
    solver.solve(0);
    time.stop();

    comm.printf("Time elapsed = %s\n", time.elapsedTime().c_str());

    viewer.write(1, comm);

    Communicator::finalize();

    return 0;
}

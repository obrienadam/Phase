#include <iostream>

#include "Input.h"
#include "CommandLine.h"
#include "ConstructGrid.h"
#include "CgnsViewer.h"
#include "Poisson.h"
#include "Time.h"

int main(int argc, char *argv[])
{
    using namespace std;

    Communicator::init(argc, argv);

    Input input;
    CommandLine(argc, argv);

    input.parseInputFile();

    shared_ptr<FiniteVolumeGrid2D> grid(constructGrid(input));
    grid->partition(input, std::make_shared<Communicator>());

    Poisson solver(input, grid);
    CgnsViewer viewer(input, solver);

    Time time;

    time.start();
    solver.solve(0);
    time.stop();

    solver.printf("Time elapsed = %s\n", time.elapsedTime().c_str());

    viewer.write(1);

    Communicator::finalize();

    return 0;
}

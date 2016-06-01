#include <iostream>

#include "Input.h"
#include "CommandLine.h"
#include "ConstructGrid.h"
#include "TecplotViewer.h"
#include "Poisson.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    Input input;
    CommandLine(argc, argv, input);

    input.parseInputFile();

    shared_ptr<FiniteVolumeGrid2D> gridPtr(constructGrid(input));

    Poisson solver(*gridPtr, input);
    TecplotViewer viewer(solver, input);

    solver.solve(0.);
    viewer.write(0.);

    return 0;
}

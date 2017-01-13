#include <iostream>

#include "Input.h"
#include "CommandLine.h"
#include "ConstructGrid.h"
#include "FractionalStep.h"
#include "EigenSparseMatrixSolver.h"
#include "Viewer.h"
#include "RunControl.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    Input input;
    CommandLine(argc, argv);

    input.parseInputFile();

    shared_ptr<FiniteVolumeGrid2D> gridPtr(constructGrid(input));
    FractionalStep solver(input, *gridPtr);

    Viewer viewer(solver, input);
    RunControl runControl;

    runControl.run(input, solver, viewer);
}


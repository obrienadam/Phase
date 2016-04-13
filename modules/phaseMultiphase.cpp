#include <iostream>

#include "Input.h"
#include "CommandLine.h"
#include "ConstructGrid.h"
#include "Multiphase.h"
#include "TecplotViewer.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    Input input;
    CommandLine(argc, argv, input);

    input.parseInputFile();

    shared_ptr<FiniteVolumeGrid2D> gridPtr(constructGrid(input));
    Multiphase solver(*gridPtr, input);
    TecplotViewer viewer(solver, input);

    Scalar maxTime = input.caseInput().get<Scalar>("Solver.maxTime");
    Scalar timeStep = input.caseInput().get<Scalar>("Solver.timeStep");
    Scalar time = 0.;

    for(time = 0.; time < maxTime; time += timeStep)
    {
        viewer.write(time);
        solver.solve(timeStep);
    }

    viewer.write(time);

    return 0;
}


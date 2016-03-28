#include <iostream>

#include "Input.h"
#include "CommandLine.h"
#include "ConstructGrid.h"
#include "Simple.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    Input input;
    CommandLine(argc, argv, input);

    input.parseInputFile();

    shared_ptr<FiniteVolumeGrid2D> gridPtr(constructGrid(input));
    Simple solver(*gridPtr, input);

    Scalar maxTime = input.caseInput().get<Scalar>("Solver.maxTime");
    Scalar timeStep = input.caseInput().get<Scalar>("Solver.timeStep");

    for(Scalar time = 0.; time < maxTime; time += timeStep)
    {
        solver.solve(timeStep);
    }

    return 0;
}


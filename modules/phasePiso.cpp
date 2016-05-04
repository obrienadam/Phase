#include <iostream>

#include "Input.h"
#include "CommandLine.h"
#include "ConstructGrid.h"
#include "Piso.h"
#include "TecplotViewer.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    Input input;
    CommandLine(argc, argv, input);

    input.parseInputFile();

    shared_ptr<FiniteVolumeGrid2D> gridPtr(constructGrid(input));
    Piso solver(*gridPtr, input);
    TecplotViewer viewer(solver, input);

    Scalar maxTime = input.caseInput().get<Scalar>("Solver.maxTime");
    Scalar timeStep = input.caseInput().get<Scalar>("Solver.timeStep");
    Scalar prevTimeStep = timeStep;

    Scalar time;

    viewer.write(time = 0.);
    solver.u.save(2);
    for(; time < maxTime; time += timeStep)
    {
        solver.solve(timeStep, prevTimeStep);
        viewer.write(time);
        printf("Simulation time: %.2lf s (%.2lf%% complete.)\n", time, time/maxTime*100);
        prevTimeStep = timeStep;
    }

    viewer.write(time);

    return 0;
}


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
    Scalar maxCo = input.caseInput().get<Scalar>("Solver.maxCo");
    Scalar time = 0.;

    Scalar maxTimeStep = input.caseInput().get<Scalar>("Solver.timeStep");
    Scalar timeStep = maxTimeStep;

    size_t fileWriteFrequency = input.caseInput().get<size_t>("System.fileWriteFrequency"), iterNo;

    for(
        time = 0., iterNo = 0;
        time < maxTime;
        time += timeStep, timeStep = std::min(solver.computeMaxTimeStep(maxCo), maxTimeStep), ++iterNo
        )
    {
        if(iterNo%fileWriteFrequency == 0)
            viewer.write(time);

        solver.solve(timeStep);
        printf("Simulation time: %.2lf s (%.2lf%% complete.)\n", time, time/maxTime*100);
    }

    viewer.write(time);

    return 0;
}


#include "RunControl.h"

void RunControl::run(const Input &input, Solver &solver, Viewer &viewer)
{
    //- Time step conditions
    Scalar maxTime = input.caseInput().get<Scalar>("Solver.maxTime");
    Scalar maxCo = input.caseInput().get<Scalar>("Solver.maxCo");

    //- Time
    Scalar time = 0.;
    Scalar maxTimeStep = input.caseInput().get<Scalar>("Solver.timeStep");
    Scalar timeStep = maxTimeStep;

    //- Write control
    size_t fileWriteFrequency = input.caseInput().get<size_t>("System.fileWriteFrequency"), iterNo;

    //- Initial conditions
    solver.setInitialConditions(input);

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
}

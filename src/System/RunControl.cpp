#include "RunControl.h"

void RunControl::run(const CommandLine &cl,
                     const Input &input,
                     SolverInterface &solver,
                     PostProcessingInterface &postProcessing)
{
    //- Max run time (for easy restart)
    auto maxWallTime =
            input.caseInput().get<Scalar>("Solver.maxWallTime", std::numeric_limits<Scalar>::infinity()) * 3600;

    //- Time step conditions
    Scalar maxTime = input.caseInput().get<Scalar>("Solver.maxTime");
    Scalar maxCo = input.caseInput().get<Scalar>("Solver.maxCo");

    //- Print the solver info
    solver.printf("%s\n", (std::string(96, '-')).c_str());
    solver.printf("%s", solver.info().c_str());
    solver.printf("%s\n", (std::string(96, '-')).c_str());

    //- Initial conditions
    solver.setInitialConditions(cl, input);
    solver.initialize();

    //- Time
    Scalar time = solver.getStartTime();
    Scalar timeStep = input.caseInput().get<Scalar>("Solver.initialTimeStep", solver.maxTimeStep());

    //- Initial output
    postProcessing.compute(0., true);

    time_.start();
    for (
            size_t iterNo = 0;
            time < maxTime && time_.elapsedSeconds(solver.comm()) < maxWallTime;
            time += timeStep, timeStep = solver.computeMaxTimeStep(maxCo, timeStep), ++iterNo
            )
    {
        solver.solve(timeStep);
        postProcessing.compute(time + timeStep, false);

        time_.stop();

        solver.printf("Time step: %.2e s\n", timeStep);
        solver.printf("Simulation time: %.2lf s (%.2lf%% complete.)\n", time + timeStep,
                      (time + timeStep) / maxTime * 100);
        solver.printf("Elapsed time: %s\n", time_.elapsedTime().c_str());
        solver.printf("Average time per iteration: %.2lf s.\n", time_.elapsedSeconds() / (iterNo + 1));
        solver.printf("%s\n", (std::string(96, '-') + "| End of iteration no " + std::to_string(iterNo + 1)).c_str());
    }
    time_.stop();

    solver.printf("%s\n", (std::string(96, '*')).c_str());
    solver.printf("Calculation complete.\n");
    solver.printf("Elapsed time: %s\n", time_.elapsedTime().c_str());
    solver.printf("Elapsed CPU time: %s\n", time_.elapsedCpuTime(solver.comm()).c_str());
    solver.printf("%s\n", (std::string(96, '*')).c_str());
}

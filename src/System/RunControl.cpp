#include "RunControl.h"

void RunControl::run(const Input &input, const Communicator &comm, Solver &solver, Viewer &viewer)
{
    //- Time step conditions
    Scalar maxTime = input.caseInput().get<Scalar>("Solver.maxTime");
    Scalar maxCo = input.caseInput().get<Scalar>("Solver.maxCo");

    //- Time
    Scalar time = 0.;
    Scalar timeStep = input.caseInput().get<Scalar>("Solver.initialTimeStep", solver.maxTimeStep());

    //- Write control
    size_t fileWriteFrequency = input.caseInput().get<size_t>("System.fileWriteFrequency"), iterNo;

    //- Print the solver info
    comm.printf("%s\n", (std::string(96, '-')).c_str());
    comm.printf("%s", solver.info().c_str());
    comm.printf("%s\n", (std::string(96, '-')).c_str());

    //- Initial conditions
    solver.setInitialConditions(input);
    solver.initialize();

    time_.start();
    for(
        time = 0., iterNo = 0;
        time < maxTime;
        time += timeStep, timeStep = solver.computeMaxTimeStep(maxCo, timeStep), ++iterNo
        )
    {
        if(iterNo%fileWriteFrequency == 0)
        {
            viewer.write(time, comm);

          //  for(VolumeIntegrator &vi: solver.volumeIntegrators())
          //      vi.integrate();

          //  for(ForceIntegrator &fi: solver.forceIntegrators())
          //      fi.integrate();

          //  viewer.write(solver.volumeIntegrators());
        }

        solver.solve(timeStep);
        time_.stop();
        comm.printf("Time step: %.2e s\n", timeStep);
        comm.printf("Simulation time: %.2lf s (%.2lf%% complete.)\n", time, time/maxTime*100);
        comm.printf("Average time per iteration: %.2lf s.\n", time_.elapsedSeconds()/(iterNo + 1));
        comm.printf("%s\n", (std::string(96, '-') + "| End of iteration no " + std::to_string(iterNo + 1)).c_str());
    }
    time_.stop();

    viewer.write(time, comm);
    comm.printf("%s\n", (std::string(96, '*')).c_str());
    comm.printf("Calculation complete.\n");
    comm.printf("Elapsed time: %s\n", time_.elapsedTime().c_str());
    comm.printf("Elapsed CPU time: %s\n", time_.elapsedCpuTime(comm).c_str());
    comm.printf("%s\n", (std::string(96, '*')).c_str());
}

#include "System/RunControl.h"

#include "Solvers/SolverFactory.h"
#include "PostProcessing/PostProcessing.h"

int main(int argc, char *argv[])
{
    using namespace std;

    Communicator::init(argc, argv);

    Input input;
    input.parseInputFile();

    std::shared_ptr<Solver> solver = SolverFactory::create(input);

    PostProcessing postProcessing(input, *solver);

    RunControl runControl;

    runControl.run(input, *solver, postProcessing);

    Communicator::finalize();
}
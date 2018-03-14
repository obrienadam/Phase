#include "SolverFactory.h"
#include "RunControl.h"

int main(int argc, char *argv[])
{
    using namespace std;

    Communicator::init(argc, argv);

    Input input;
    input.parseInputFile();

    std::shared_ptr<Solver> solver = SolverFactory::create(input);
    CgnsViewer viewer(input, *solver);
    RunControl runControl;

    runControl.run(input, *solver, viewer);

    Communicator::finalize();
}
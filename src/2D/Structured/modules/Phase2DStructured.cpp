#include "System/CommandLine.h"
#include "System/RunControl.h"

#include "StructuredGrid2D/StructuredGrid2D.h"
#include "Solvers/SolverFactory.h"
#include "PostProcessing/PostProcessing.h"

int main(int argc, char *argv[])
{
    Communicator::init(argc, argv);

    CommandLine cl;

    Input input;
    input.parseInputFile();

    auto grid = std::make_shared<StructuredGrid2D>(input);

    auto solver = SolverFactory::create(input, grid);

    auto postProcessing = PostProcessing(input, *solver);

    RunControl runControl;

    //runControl.run(cl, input, *solver, postProcessing);

    Communicator::finalize();

    return 0;
}

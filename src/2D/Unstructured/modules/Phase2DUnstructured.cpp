#include "System/CommandLine.h"
#include "System/RunControl.h"

#include "FiniteVolumeGrid2D/FiniteVolumeGrid2DFactory.h"
#include "Solvers/SolverFactory.h"
#include "PostProcessing/PostProcessing.h"

int main(int argc, char *argv[])
{
    using namespace std;

    Communicator::init(argc, argv);

    CommandLine cl;

    cl.addSwitch("restart,r", "restart the solution from the latest time point");
    cl.addSwitch("use-partitioned-grid,g", "use the pre-partitioned grid");

    cl.parseArguments(argc, argv);

    Input input;
    input.parseInputFile();

    std::shared_ptr<FiniteVolumeGrid2D> grid = FiniteVolumeGrid2DFactory::create(cl, input);

    std::shared_ptr<Solver> solver = SolverFactory::create(input, grid);

    PostProcessing postProcessing(input, *solver, std::weak_ptr<const ImmersedBoundary>());

    RunControl runControl;

    runControl.run(cl, input, *solver, postProcessing);

    Communicator::finalize();
}

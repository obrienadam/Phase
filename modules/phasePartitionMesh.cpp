#include <map>
#include <metis.h>

#include "CommandLine.h"
#include "Input.h"
#include "CgnsUnstructuredGrid.h"
#include "Communicator.h"
#include "Poisson.h"
#include "Viewer.h"
#include "Exception.h"

int main(int argc, char *argv[])
{
    using namespace std;

    Communicator::init(argc, argv);
    Communicator comm(MPI_COMM_WORLD);

    Input input;
    input.parseInputFile();
    CgnsUnstructuredGrid grid(input);
    grid.partition(comm);
    grid.save("test.cgns", comm);

    Poisson solver(input, comm, grid);

    for(Scalar &val: solver.phi)
        val = comm.rank();

    grid.sendMessages(comm, solver.phi);

    if(comm.isMainProc())
    {
        Viewer viewer(input, solver);
        viewer.write(0);
    }

    Communicator::finalize();

    return 0;
}

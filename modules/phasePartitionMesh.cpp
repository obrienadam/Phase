#include <map>
#include <metis.h>

#include "CommandLine.h"
#include "Input.h"
#include "CgnsUnstructuredGrid.h"
#include "Communicator.h"
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

    Communicator::finalize();
    return 0;
}

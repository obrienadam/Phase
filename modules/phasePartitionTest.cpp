#include <memory.h>

#include "Communicator.h"
#include "CgnsUnstructuredGrid.h"

int main(int argc, char *argv[])
{
    Communicator::init();

    Communicator comm(MPI_COMM_WORLD);
    Input input;

    input.parseInputFile();
    CgnsUnstructuredGrid grid(input);

    comm.printf("%s\n", "Hello world!");
    comm.partitionGrid(grid);

    Communicator::finalize();
}

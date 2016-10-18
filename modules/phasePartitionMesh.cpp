#include <map>
#include <metis.h>

#include "CommandLine.h"
#include "Input.h"
#include "CgnsUnstructuredGrid.h"
#include "Exception.h"

int main(int argc, const char *argv[])
{
    using namespace std;

    CommandLine cl;

    cl.setOptions(map<string, string>{
                      {"--np", "Number of mesh partitions"},
                  });

    cl.parseArguments(argc, argv);

    Input input;
    input.parseInputFile();

    CgnsUnstructuredGrid grid(input);
    CgnsUnstructuredGrid partitionedGrid();

    int nElems = grid.nCells();
    int nNodes = grid.nNodes();

    vector<idx_t> elems;
    vector<idx_t> elemInds;

    elems.reserve(4*nElems);
    elemInds.reserve(nElems + 1);
    elemInds.push_back(0); // for the first element

    for(const Cell &cell: grid.cells())
    {
        const Size nVerts = cell.nodes().size();
        elemInds.push_back(elemInds.back() + nVerts);

        for(const Node &node: cell.nodes())
            elems.push_back(node.id());
    }

    idx_t nCommon = 2;
    idx_t nPartitions = std::stoi(cl.getOption("--np"));
    idx_t objVal;

    vector<idx_t> elemPart(grid.nCells());
    vector<idx_t> nodePart(grid.nNodes());

    printf("Partitioning mesh into %d subdomains...\n", nPartitions);
    int status = METIS_PartMeshDual(&nElems, &nNodes,
                                    elemInds.data(), elems.data(),
                                    NULL, NULL,
                                    &nCommon, &nPartitions,
                                    NULL, NULL, &objVal,
                                    elemPart.data(), nodePart.data());
    if(status == METIS_OK)
        printf("Sucessfully partitioned mesh.\n");
    else
        throw Exception("phasePartitionMesh", "phasePartitionMesh", "an error occurred during partitioning.");

    return 0;
}

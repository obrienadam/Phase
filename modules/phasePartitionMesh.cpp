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

    vector<idx_t> cellPartition(grid.nCells());
    vector<idx_t> nodePartition(grid.nNodes());

    printf("Partitioning mesh into %d subdomains...\n", nPartitions);
    int status = METIS_PartMeshDual(&nElems, &nNodes,
                                    elemInds.data(), elems.data(),
                                    NULL, NULL,
                                    &nCommon, &nPartitions,
                                    NULL, NULL, &objVal,
                                    cellPartition.data(), nodePartition.data());
    if(status == METIS_OK)
        printf("Sucessfully partitioned mesh.\n");
    else
        throw Exception("phasePartitionMesh", "phasePartitionMesh", "an error occurred during partitioning.");

    //- Construct a cell list for each partition
    vector< vector<Label> > partitionCellLists(nPartitions);

    for(const Cell& cell: grid.cells())
        partitionCellLists[cellPartition[cell.id()]].push_back(cell.id());

    //- Create the new file for the partitioned grid
    CgnsUnstructuredGrid partitionedGrid;
    partitionedGrid.newFile("partitionedMesh.cgns", input.caseInput().get<string>("CaseName"));

    //- Loop over all cell lists. Construct the nodes and cells for each partition
    for(int partitionNo = 0; partitionNo < nPartitions; ++partitionNo)
    {
        const vector<Ref<const Cell>> cellList = grid.getCells(partitionCellLists[partitionNo]);
        map<Label, Label> nodeIdMap;
        vector<Point2D> nodes;
        vector<cgsize_t> cells;
        cells.reserve(5*cellList.size());

        //- Construct node list and node map, as well as cells
        Label locId = 0;
        for(const Cell& cell: cellList)
        {
            cells.push_back(cell.nodes().size() == 3 ? TRI_3 : QUAD_4);

            for(const Node& node: cell.nodes())
            {
                auto it = nodeIdMap.find(node.id());

                if(it == nodeIdMap.end())
                {
                    nodeIdMap[node.id()] = locId++;
                    nodes.push_back(node);
                }

                cells.push_back(nodeIdMap[node.id()] + 1);
            }
        }

        int zoneNo = partitionedGrid.addZone("partition_" + to_string(partitionNo), nodes.size(), cellList.size());
        partitionedGrid.addNodes(zoneNo, nodes);
        partitionedGrid.addMixedCells(zoneNo, cellList.size(), cells);
    }

    partitionedGrid.closeFile();

    return 0;
}

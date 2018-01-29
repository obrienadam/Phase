#include <iostream>
#include <unordered_set>

#include <boost/filesystem.hpp>
#include <metis.h>
#include <cgnslib.h>

#include "Input.h"
#include "ConstructGrid.h"

int main(int argc, char *argv[])
{
    using namespace std;

    Input input;
    input.parseInputFile();

    auto grid = constructGrid(input);

    cout << "Successfully loaded grid.\n";

    Scalar minBufferWidth = input.caseInput().get<Scalar>("Grid.minBufferWidth", 0.);

    idx_t nPartitions = 32;
    idx_t nElems = grid->nCells();
    pair<vector<int>, vector<int>> mesh = grid->nodeElementConnectivity();
    idx_t nNodes = grid->nNodes();
    idx_t nCommon = 2; //- face connectivity weighting only
    idx_t objVal;
    vector<idx_t> nodePartition(nNodes);
    vector<idx_t> cellPartition(nElems);

    int status = METIS_PartMeshDual(&nElems, &nNodes,
                                    mesh.first.data(), mesh.second.data(),
                                    NULL, NULL,
                                    &nCommon, &nPartitions,
                                    NULL, NULL, &objVal,
                                    cellPartition.data(), nodePartition.data());
    if (status == METIS_OK)
        cout << "Successfully computed partitioning. Objective val = " << objVal << "\n";
    else
        throw Exception("", "", "an error occurred during partitioning.");

    vector<unordered_set<int>> cellsOnProc(nPartitions);
    vector<vector<unordered_set<int>>> bufferCells(nPartitions, vector<unordered_set<int>>(nPartitions));

    for(const Cell& cell: grid->cells())
    {
        cellsOnProc[cellPartition[cell.id()]].insert(cell.id());

        for(const CellLink& nb: cell.cellLinks())
            if(cellPartition[cell.id()] != cellPartition[nb.cell().id()])
            {
                cellsOnProc[cellPartition[nb.cell().id()]].insert(nb.cell().id());
                bufferCells[cellPartition[cell.id()]][cellPartition[nb.cell().id()]].insert(nb.cell().id());
            }
    }

    unordered_map<string, vector<int>> patches = grid->patchToNodeMap();

    grid = nullptr;

    for(int proc = 0; proc < nPartitions; ++proc)
    {
        const auto &cells = cellsOnProc[proc];
        const auto &buffCells = bufferCells[proc];

        vector<double> coordsX, coordsY;
        unordered_map<int, int> globalToLocalNode;

        for(int cell: cells)
        {
            for(int node = mesh.first[cell]; node < mesh.first[cell + 1]; ++node)
            {

            }
        }


        int fid;
        cg_open(("solution/Proc" + std::to_string(proc) + "Grid.cgns").c_str(), CG_MODE_WRITE, &fid);

        int cid;
        cg_coord_write(fid, 1, 1, CGNS_ENUMV(RealDouble), "CoordinateX", coordsX.data(), &cid);
        cg_coord_write(fid, 1, 1, CGNS_ENUMV(RealDouble), "CoordinateY", coordsY.data(), &cid);



        cg_close(fid);
    }

    return 0;
}
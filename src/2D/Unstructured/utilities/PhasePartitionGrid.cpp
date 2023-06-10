#include <iostream>

#include <boost/filesystem.hpp>
#include <metis.h>

#include "System/CgnsFile.h"
#include "System/CommandLine.h"
#include "System/Input.h"

#include "FiniteVolumeGrid2D/FiniteVolumeGrid2DFactory.h"

int main(int argc, char *argv[]) {
  namespace po = boost::program_options;

  Communicator::init(argc, argv);

  CommandLine cl;

  cl.addOptions()("num-partitions,n", po::value<int>()->required(),
                  "Number of partitions to generate")(
      "min-buffer-width,m", po::value<double>()->required(),
      "Minimum cell buffer width");

  cl.parseArguments(argc, argv);

  idx_t numPartitions = cl.get<int>("num-partitions");
  double minBufferWidth = cl.get<double>("min-buffer-width");

  std::cout << "Attempting to partition grid into " << numPartitions
            << " partitions with a minimum cell buffer width of "
            << minBufferWidth << " units.\n";

  Input input;
  input.parseInputFile();

  std::cout << "Generating grid...\n";

  auto grid = FiniteVolumeGrid2DFactory::create(input);

  std::cout << "Computing partitioning...\n";

  idx_t ncon = 1, nvtxs = grid->nCells(), nparts = numPartitions, objval;

  auto graph = grid->connectivityGraph();

  std::vector<idx_t> cellPartition(grid->nCells());

  int status = METIS_PartGraphRecursive(
      &nvtxs, &ncon, graph.first.data(), graph.second.data(), NULL, NULL, NULL,
      &nparts, NULL, NULL, NULL, &objval, cellPartition.data());
  if (status == METIS_OK)
    std::cout << "Sucessfully computed partitioning. Objective value = "
              << objval << "\n";
  else
    throw Exception("", "PhasePartitionGrid",
                    "an error occurred during partitioning.");

  //- Create cgns files for each partition
  for (int proc = 0; proc < numPartitions; ++proc) {
    std::cout << "Creating grid file for proc " << proc << ".\n";

    std::vector<Label> localCells;
    std::vector<Label> partitionBoundaryCells;
    std::vector<int> owningProc;

    for (const Cell &cell : grid->cells())
      if (cellPartition[cell.id()] == proc) {
        localCells.push_back(cell.id());
        owningProc.push_back(proc);

        for (const CellLink &nb : cell.cellLinks())
          if (cellPartition[nb.cell().id()] != proc) {
            partitionBoundaryCells.push_back(cell.id());
            break;
          }
      }

    std::unordered_set<Label> tmp;

    for (Label id : partitionBoundaryCells) {
      const Cell &cell = grid->cells()[id];

      for (const CellLink &nb : cell.cellLinks())
        if (cellPartition[nb.cell().id()] != proc &&
            tmp.insert(nb.cell().id()).second) {
          localCells.push_back(nb.cell().id());
          owningProc.push_back(cellPartition[nb.cell().id()]);
        }

      //- Expands the buffer region
      for (const Cell &kcell : grid->localCells().itemsCoveredBy(
               Circle(cell.centroid(), minBufferWidth)))
        if (cellPartition[kcell.id()] != proc &&
            tmp.insert(kcell.id()).second) {
          localCells.push_back(kcell.id());
          owningProc.push_back(cellPartition[kcell.id()]);
        }
    }

    std::unordered_map<Label, Label> globalToLocalNodeId;
    std::vector<Point2D> localNodes;
    std::vector<int> eptr(1, 0), eind;

    for (Label id : localCells) {
      const Cell &cell = grid->cells()[id];
      eptr.push_back(eptr.back() + cell.nodes().size());

      for (const Node &node : cell.nodes()) {
        auto insert = globalToLocalNodeId.insert(
            std::make_pair(node.id(), localNodes.size()));

        if (insert.second)
          localNodes.push_back(node);

        eind.push_back(insert.first->second + 1);
      }
    }

    std::unordered_map<std::string, std::vector<int>> patches;

    for (const FaceGroup &patch : grid->patches())
      for (const Face &face : patch) {
        auto itr1 = globalToLocalNodeId.find(face.lNode().id());
        auto itr2 = globalToLocalNodeId.find(face.rNode().id());

        if (itr1 != globalToLocalNodeId.end() &&
            itr2 != globalToLocalNodeId.end()) {
          patches[patch.name()].push_back(itr1->second + 1);
          patches[patch.name()].push_back(itr2->second + 1);
        }
      }

    //- Make sure the output directory is available
    boost::filesystem::path path = "./solution/Proc" + std::to_string(proc);
    boost::filesystem::create_directories(path);

    CgnsFile file((path / "Grid.cgns").string(), CgnsFile::WRITE);

    int bid =
        file.createBase(input.caseInput().get<std::string>("CaseName"), 2, 2);
    int zid = file.createUnstructuredZone(bid, "Zone", localNodes.size(),
                                          localCells.size());

    file.writeCoordinates(bid, zid, localNodes);

    int start = 1;
    int end = localCells.size();

    file.writeMixedElementSection(bid, zid, "Cells", start, end, eptr, eind);

    for (const auto &entry : patches) {
      start = end + 1;
      end = start + entry.second.size() / 2 - 1;

      file.writeBarElementSection(bid, zid, entry.first, start, end,
                                  entry.second);
      file.writeBoCo(bid, zid, entry.first, start, end);
    }

    int sid = file.writeSolution(bid, zid, "Info");

    file.writeField(bid, zid, sid, "GlobalID", localCells);
    file.writeField(bid, zid, sid, "ProcNo", owningProc);

    file.close();
  }

  Communicator::finalize();
}
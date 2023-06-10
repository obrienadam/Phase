#include <fstream>

#include "CgnsUnstructuredGrid.h"
#include "FiniteVolumeGrid2DFactory.h"
#include "StructuredRectilinearGrid.h"

std::shared_ptr<FiniteVolumeGrid2D>
FiniteVolumeGrid2DFactory::create(GridType type, const Input &input) {
  std::shared_ptr<FiniteVolumeGrid2D> grid = nullptr;

  switch (type) {
  case RECTILINEAR:
    grid = std::make_shared<StructuredRectilinearGrid>(input);
    break;
  case CGNS:
    grid = std::make_shared<CgnsUnstructuredGrid>(input);
    break;
  case COORDS: {
    std::vector<Scalar> xcoords, ycoords;
    std::ifstream fin;

    fin.open(input.caseInput().get<std::string>("Grid.xCoordFile"));

    if (!fin.is_open())
      throw Exception("FiniteVolumeGrid2DFactory", "create",
                      "missing x-coordinate file.");

    std::string line;
    while (!fin.eof()) {
      std::getline(fin, line);
      xcoords.emplace_back(std::stod(line));
    }
    fin.close();
    fin.open(input.caseInput().get<std::string>("Grid.yCoordFile"));

    if (!fin.is_open())
      throw Exception("FiniteVolumeGrid2DFactory", "create",
                      "missing y-coordinate file.");

    while (!fin.eof()) {
      std::getline(fin, line);
      ycoords.emplace_back(std::stod(line));
    }
    fin.close();

    grid = std::make_shared<StructuredRectilinearGrid>(xcoords, ycoords);
  } break;
  case LOAD:
    auto grid = std::make_shared<CgnsUnstructuredGrid>();
    grid->load("./solution/Proc" + std::to_string(grid->comm().rank()) +
                   "/Grid.cgns",
               Vector2D(0., 0.));
    grid->readPartitionData("./solution/Proc" +
                            std::to_string(grid->comm().rank()) + "/Grid.cgns");
    return grid;
  }

  grid->partition(input);

  return grid;
}

std::shared_ptr<FiniteVolumeGrid2D>
FiniteVolumeGrid2DFactory::create(std::string type, const Input &input) {
  std::transform(type.begin(), type.end(), type.begin(),
                 [](unsigned char c) { return std::tolower(c); });

  if (type == "rectilinear")
    return create(RECTILINEAR, input);
  else if (type == "cgns")
    return create(CGNS, input);
  else if (type == "coords")
    return create(COORDS, input);
  else if (type == "load")
    return create(LOAD, input);

  throw Exception("FiniteVolumeGrid2DFactory", "create",
                  "grid \"" + type + "\" is not a valid grid type.");
}

std::shared_ptr<FiniteVolumeGrid2D>
FiniteVolumeGrid2DFactory::create(const Input &input) {
  return create(input.caseInput().get<std::string>("Grid.type"), input);
}

std::shared_ptr<FiniteVolumeGrid2D>
FiniteVolumeGrid2DFactory::create(const CommandLine &cl, const Input &input) {
  if (cl.get<bool>("use-partitioned-grid") || cl.get<bool>("restart"))
    return create(LOAD, input);
  else
    return create(input);
}

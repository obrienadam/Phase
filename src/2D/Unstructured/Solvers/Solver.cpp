#include <math.h>
#include <regex>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include "System/CgnsFile.h"

#include "Solver.h"

Solver::Solver(const Input &input,
               const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
    : grid_(grid) {
  //- Set simulation time options
  maxTimeStep_ = input.caseInput().get<Scalar>("Solver.timeStep");
  startTime_ = 0.;

  //- Index map
  scalarIndexMap_ = std::make_shared<IndexMap>(*grid_, 1);
  vectorIndexMap_ = std::make_shared<IndexMap>(*grid_, 2);
}

int Solver::printf(const char *format, ...) const {
  int n = 0;

  if (comm().isMainProc()) {
    va_list argsPtr;
    va_start(argsPtr, format);
    n = vfprintf(stdout, format, argsPtr);
    va_end(argsPtr);
  }

  return n;
}

Scalar Solver::getStartTime() const { return startTime_; }

std::shared_ptr<FiniteVolumeField<int>>
Solver::integerField(const std::string &name) const {
  auto it = integerFields_.find(name);
  return it != integerFields_.end() ? it->second : nullptr;
}

std::shared_ptr<ScalarFiniteVolumeField>
Solver::scalarField(const std::string &name) const {
  auto it = scalarFields_.find(name);
  return it != scalarFields_.end() ? it->second : nullptr;
}

std::shared_ptr<VectorFiniteVolumeField>
Solver::vectorField(const std::string &name) const {
  auto it = vectorFields_.find(name);
  return it != vectorFields_.end() ? it->second : nullptr;
}

template <>
std::shared_ptr<FiniteVolumeField<int>>
Solver::addField(const std::string &name,
                 const std::shared_ptr<CellGroup> &cells) {
  auto insert = integerFields_.insert(std::make_pair(
      name, std::make_shared<FiniteVolumeField<int>>(
                grid_, name, 0, true, false, cells, scalarIndexMap_)));

  if (!insert.second)
    throw Exception("Solver", "addField",
                    "field \"" + name + "\" already exists.");

  return insert.first->second;
}

template <>
std::shared_ptr<FiniteVolumeField<Scalar>>
Solver::addField(const std::string &name,
                 const std::shared_ptr<CellGroup> &cells) {
  auto insert = scalarFields_.insert(std::make_pair(
      name, std::make_shared<ScalarFiniteVolumeField>(
                grid_, name, 0, true, false, cells, scalarIndexMap_)));

  if (!insert.second)
    throw Exception("Solver", "addField",
                    "field \"" + name + "\" already exists.");

  return insert.first->second;
}

template <>
std::shared_ptr<FiniteVolumeField<Vector2D>>
Solver::addField(const std::string &name,
                 const std::shared_ptr<CellGroup> &cells) {
  auto insert = vectorFields_.insert(std::make_pair(
      name, std::make_shared<VectorFiniteVolumeField>(
                grid_, name, Vector2D(), true, false, cells, vectorIndexMap_)));

  if (!insert.second)
    throw Exception("Solver", "addField",
                    "field \"" + name + "\" already exists.");

  return insert.first->second;
}

template <>
std::shared_ptr<FiniteVolumeField<Scalar>>
Solver::addField(const Input &input, const std::string &name,
                 const std::shared_ptr<CellGroup> &cells) {
  auto insert = scalarFields_.insert(std::make_pair(
      name, std::make_shared<ScalarFiniteVolumeField>(
                input, grid_, name, 0., true, false, cells, scalarIndexMap_)));

  if (!insert.second)
    throw Exception("Solver", "addField",
                    "field \"" + name + "\" already exists.");

  return insert.first->second;
}

template <>
std::shared_ptr<FiniteVolumeField<Vector2D>>
Solver::addField(const Input &input, const std::string &name,
                 const std::shared_ptr<CellGroup> &cells) {
  auto insert = vectorFields_.insert(std::make_pair(
      name, std::make_shared<VectorFiniteVolumeField>(input, grid_, name,
                                                      Vector2D(), true, false,
                                                      cells, vectorIndexMap_)));

  if (!insert.second)
    throw Exception("Solver", "addField",
                    "field \"" + name + "\" already exists.");

  return insert.first->second;
}

template <>
std::shared_ptr<FiniteVolumeField<int>>
Solver::addField(const std::shared_ptr<FiniteVolumeField<int>> &field) {
  auto insert = integerFields_.insert(std::make_pair(field->name(), field));

  if (!insert.second)
    throw Exception("Solver", "addField",
                    "field \"" + field->name() + "\" already exists.");

  // insert.first->second->setIndexMap(scalarIndexMap_);

  return insert.first->second;
}

template <>
std::shared_ptr<FiniteVolumeField<Scalar>>
Solver::addField(const std::shared_ptr<FiniteVolumeField<Scalar>> &field) {
  auto insert = scalarFields_.insert(std::make_pair(field->name(), field));

  if (!insert.second)
    throw Exception("Solver", "addField",
                    "field \"" + field->name() + "\" already exists.");

  insert.first->second->setIndexMap(scalarIndexMap_);

  return insert.first->second;
}

template <>
std::shared_ptr<FiniteVolumeField<Vector2D>>
Solver::addField(const std::shared_ptr<FiniteVolumeField<Vector2D>> &field) {
  auto insert = vectorFields_.insert(std::make_pair(field->name(), field));

  if (!insert.second)
    throw Exception("Solver", "addField",
                    "field \"" + field->name() + "\" already exists.");

  insert.first->second->setIndexMap(vectorIndexMap_);

  return insert.first->second;
}

template <>
std::shared_ptr<FiniteVolumeField<Tensor2D>>
Solver::addField(const std::shared_ptr<FiniteVolumeField<Tensor2D>> &field) {
  auto insert = tensorFields_.insert(std::make_pair(field->name(), field));

  if (!insert.second)
    throw Exception("Solver", "addField",
                    "field \"" + field->name() + "\" already exists.");

  // insert.first->second->setIndexMap(indexMap_);

  return insert.first->second;
}

void Solver::setInitialConditions(const Input &input) {
  using namespace std;
  using namespace boost::property_tree;

  for (const auto &child :
       input.initialConditionInput().get_child("InitialConditions")) {
    auto scalarFieldIt = scalarFields_.find(child.first);

    if (scalarFieldIt != scalarFields_.end()) {
      ScalarFiniteVolumeField &field = *scalarFieldIt->second;

      for (const auto &ic : child.second) {
        const auto &icTree = ic.second;
        std::string type = icTree.get<string>("type");

        if (type == "circle") {
          Circle circle = Circle(Vector2D(icTree.get<string>("center")),
                                 icTree.get<Scalar>("radius"));
          setCircle(circle, icTree.get<Scalar>("value"), field);
        } else if (type == "circleSector") {
          Circle circle = Circle(Vector2D(icTree.get<string>("center")),
                                 icTree.get<Scalar>("radius"));
          Scalar thetaMin = icTree.get<Scalar>("thetaMin");
          Scalar thetaMax = icTree.get<Scalar>("thetaMax");
          Scalar innerValue = icTree.get<Scalar>("value");
          setCircleSector(circle, thetaMin, thetaMax, innerValue, field);
        } else if (type == "box") {
          Point2D center = Point2D(icTree.get<string>("center"));
          Scalar w = icTree.get<Scalar>("width") / 2;
          Scalar h = icTree.get<Scalar>("height") / 2;

          Polygon pgn = {Point2D(center.x - w, center.y - h),
                         Point2D(center.x + w, center.y - h),
                         Point2D(center.x + w, center.y + h),
                         Point2D(center.x - w, center.y + h)};

          setBox(pgn, icTree.get<Scalar>("value"), field);
        } else if (type == "uniform")
          field.fillInterior(icTree.get<Scalar>("value"));
        else if (type == "rotating") {
          setRotating(icTree.get<std::string>("function"),
                      icTree.get<Scalar>("amplitude"),
                      Vector2D(icTree.get<std::string>("center")), field);
        } else
          throw Exception("Input", "setInitialConditions",
                          "invalid initial condition type \"" + type + "\".");

        printf("Set initial condition \"%s\" of type %s on field \"%s\".\n",
               ic.first.c_str(), type.c_str(), field.name().c_str());
      }

      continue;
    }

    auto vectorFieldIt = vectorFields_.find(child.first);

    if (vectorFieldIt != vectorFields_.end()) {
      VectorFiniteVolumeField &field = *vectorFieldIt->second;

      for (const auto &ic : child.second) {
        const auto &icTree = ic.second;
        std::string type = icTree.get<string>("type");

        if (type == "circle") {
          Circle circle = Circle(Vector2D(icTree.get<string>("center")),
                                 icTree.get<Scalar>("radius"));
          setCircle(circle, Vector2D(icTree.get<string>("value")), field);
        } else if (type == "square") {
          Point2D center = Point2D(icTree.get<string>("center"));
          Scalar w = icTree.get<Scalar>("width") / 2;
          Scalar h = icTree.get<Scalar>("height") / 2;

          Polygon pgn = {Point2D(center.x - w, center.y - h),
                         Point2D(center.x + w, center.y - h),
                         Point2D(center.x + w, center.y + h),
                         Point2D(center.x - w, center.y + h)};

          setBox(pgn, Vector2D(icTree.get<string>("value")), field);
        } else if (type == "uniform")
          field.fillInterior(Vector2D(icTree.get<string>("value")));
        else if (type == "rotating") {
          setRotating(icTree.get<std::string>("xFunction"),
                      icTree.get<std::string>("yFunction"),
                      Vector2D(icTree.get<std::string>("amplitude")),
                      Vector2D(icTree.get<std::string>("center")), field);
        } else
          throw Exception("Input", "setInitialConditions",
                          "invalid initial condition type \"" + type + "\".");

        printf("Set initial condition \"%s\" of type %s on field \"%s\".\n",
               ic.first.c_str(), type.c_str(), field.name().c_str());
      }
    }
  }
}

void Solver::setInitialConditions(const CommandLine &cl, const Input &input) {
  if (cl.get<bool>("restart"))
    restartSolution(input);
  else
    setInitialConditions(input);
}

//- Protected methods

void Solver::setCircle(const Circle &circle, Scalar innerValue,
                       ScalarFiniteVolumeField &field) {
  const Polygon pgn = circle.polygonize(1000);

  for (const Cell &cell : field.grid()->localCells()) {
    Scalar area = 0.;

    for (const Polygon &xc : intersection(pgn, cell.shape()))
      area += xc.area();

    field(cell) = innerValue * area / cell.shape().area();
  }

  grid_->sendMessages(field);
  field.interpolateFaces();
}

void Solver::setCircle(const Circle &circle, const Vector2D &innerValue,
                       VectorFiniteVolumeField &field) {
  const Polygon pgn = circle.polygonize(1000);

  for (const Cell &cell : field.grid()->localCells()) {
    Scalar area = 0.;

    for (const Polygon &xc : intersection(pgn, cell.shape()))
      area += xc.area();

    field(cell) = innerValue * area / cell.shape().area();
  }

  grid_->sendMessages(field);
  field.interpolateFaces();
}

void Solver::setCircleSector(const Circle &circle, Scalar thetaMin,
                             Scalar thetaMax, Scalar innerValue,
                             ScalarFiniteVolumeField &field) {
  thetaMin *= M_PI / 180.;
  thetaMax *= M_PI / 180.;

  while (thetaMin > 2. * M_PI)
    thetaMin -= 2. * M_PI;

  while (thetaMin < 0.)
    thetaMin += 2. * M_PI;

  while (thetaMax > 2. * M_PI)
    thetaMax -= 2. * M_PI;

  while (thetaMax < 0.)
    thetaMax += 2. * M_PI;

  std::vector<Point2D> vtx;
  const int nPts = 1000;
  const Scalar dTheta = (thetaMax - thetaMin) / (nPts - 1);

  for (int i = 0; i < nPts; ++i) {
    const Scalar theta = thetaMin + i * dTheta;
    const Vector2D rVec(circle.radius() * std::cos(theta),
                        circle.radius() * std::sin(theta));

    vtx.push_back(circle.centroid() + rVec);
  }

  Polygon pgn(vtx.begin(), vtx.end());

  for (const Cell &cell : field.grid()->localCells()) {
    Scalar area = 0.;

    for (const Polygon &xc : intersection(pgn, cell.shape()))
      area += xc.area();

    field(cell) = innerValue * area / cell.shape().area();
  }

  grid_->sendMessages(field);
  field.interpolateFaces();
}

void Solver::setBox(const Polygon &box, Scalar innerValue,
                    ScalarFiniteVolumeField &field) {
  for (const Cell &cell : field.grid()->localCells()) {
    Scalar area = 0.;

    for (const Polygon &pgn : intersection(box, cell.shape()))
      area += pgn.area();

    field(cell) = innerValue * area / cell.shape().area();
  }

  grid_->sendMessages(field);
  field.interpolateFaces();
}

void Solver::setBox(const Polygon &box, const Vector2D &innerValue,
                    VectorFiniteVolumeField &field) {
  for (const Cell &cell : field.grid()->localCells()) {
    Scalar area = 0.;

    for (const Polygon &pgn : intersection(box, cell.shape()))
      area += pgn.area();

    field(cell) = innerValue * area / cell.shape().area();
  }

  grid_->sendMessages(field);
  field.interpolateFaces();
}

void Solver::setRotating(const std::string &function, Scalar amplitude,
                         const Vector2D &center,
                         ScalarFiniteVolumeField &field) {
  std::function<Scalar(Scalar)> func;

  if (function == "sin")
    func = [](Scalar x) { return sin(x); };
  else if (function == "cos")
    func = [](Scalar x) { return cos(x); };
  else
    throw Exception("Input", "setRotating", "invalid rotation function.");

  for (const Cell &cell : field.grid()->cells()) {
    Vector2D rVec = cell.centroid() - center;
    Scalar theta = atan2(rVec.y, rVec.x);

    field[cell.id()] = amplitude * func(theta);
  }

  for (const Face &face : field.grid()->interiorFaces()) {
    Vector2D rVec = face.centroid() - center;
    Scalar theta = atan2(rVec.y, rVec.x);

    field.faces()[face.id()] = amplitude * func(theta);
  }
}

void Solver::setRotating(const std::string &xFunction,
                         const std::string &yFunction,
                         const Vector2D &amplitude, const Vector2D &center,
                         VectorFiniteVolumeField &field) {
  std::function<Scalar(Scalar)> xFunc, yFunc;

  if (xFunction == "sin")
    xFunc = [](Scalar x) { return sin(x); };
  else if (xFunction == "cos")
    xFunc = [](Scalar x) { return cos(x); };
  else
    throw Exception("Input", "setRotating", "invalid x rotation function.");

  if (yFunction == "sin")
    yFunc = [](Scalar x) { return sin(x); };
  else if (yFunction == "cos")
    yFunc = [](Scalar x) { return cos(x); };
  else
    throw Exception("Input", "setRotating", "invalid y rotation function.");

  for (const Cell &cell : field.grid()->cells()) {
    Vector2D rVec = cell.centroid() - center;
    Scalar theta = atan2(rVec.y, rVec.x);

    field[cell.id()].x = amplitude.x * xFunc(theta);
    field[cell.id()].y = amplitude.y * yFunc(theta);
  }

  for (const Face &face : field.grid()->interiorFaces()) {
    Vector2D rVec = face.centroid() - center;
    Scalar theta = atan2(rVec.y, rVec.x);

    field.faces()[face.id()].x = amplitude.x * xFunc(theta);
    field.faces()[face.id()].y = amplitude.y * yFunc(theta);
  }
}

void Solver::restartSolution(const Input &input) {
  using namespace std;
  string viewerType = input.postProcessingInput().get<string>(
      "PostProcessing.viewerType", "cgns");

  if (viewerType == "cgns")
    restartFromCgnsViewer(input);
  else if (viewerType == "compactCgns")
    restartFromCompactCgnsViewer(input);
  else
    throw Exception("Solver", "restartSolution",
                    "Unrecognized viewer type \"" + viewerType + "\".");
}

void Solver::restartFromCgnsViewer(const Input &input) {
  using namespace std;
  using namespace boost::filesystem;

  regex re("[0-9]+\\.[0-9]+");
  path path;

  for (directory_iterator end, dir("./solution"); dir != end; ++dir)
    if (regex_match(dir->path().filename().string(), re)) {
      std::smatch match;
      std::regex_search(dir->path().filename().string(), match, re);
      Scalar time = std::stod(match.str());

      if (time >= startTime_) {
        path = dir->path();
        startTime_ = time;
      }
    }

  path /= ("Proc" + std::to_string(grid_->comm().rank())) / "Solution.cgns";

  //- Check to make sure the restart can be done
  if (!exists(path))
    throw Exception("Solver", "restartSolution",
                    "no file \"" + path.string() + "\" needed for restart.");

  startTime_ = grid_->comm().broadcast(grid_->comm().mainProcNo(), startTime_);
  grid_->comm().printf("Restarting solution at t = %lf...\n", startTime_);

  CgnsFile file(path.string(), CgnsFile::READ);
  readLatestCgnsFlowSolution(file);
  file.close();
}

void Solver::restartFromCompactCgnsViewer(const Input &input) {
  using namespace std;
  using namespace boost::filesystem;

  int proc = grid_->comm().rank();
  path path = ("./solution") / ("proc" + to_string(proc) + ".cgns");
  CgnsFile file(path.string(), CgnsFile::READ);
  int sid = readLatestCgnsFlowSolution(file);

  for (const auto &d : file.readDescriptorNodes(1, 1, sid))
    if (d.first == "SolutionTime") {
      startTime_ = stod(d.second);
      break;
    }

  startTime_ = grid_->comm().broadcast(grid_->comm().mainProcNo(), startTime_);
  grid_->comm().printf("Restarting solution at t = %lf...\n", startTime_);

  file.close();
}

int Solver::readLatestCgnsFlowSolution(const CgnsFile &file) {
  using namespace std;
  CgnsFile::Solution soln = file.readLastFlowSolution(1, 1);

  if (soln.id < 1)
    return soln.id;

  for (const auto &entry : scalarFields_) {
    auto field =
        file.readField<Scalar>(1, 1, soln.id, 1, grid_->nCells(), entry.first);

    if (field.data.size() == entry.second->size())
      copy(field.data.begin(), field.data.end(), entry.second->begin());

    grid_->sendMessages(*entry.second);
    entry.second->interpolateFaces();
  }

  for (const auto &entry : vectorFields_) {
    auto fieldX = file.readField<Scalar>(1, 1, soln.id, 1, grid_->nCells(),
                                         entry.first + "X");
    auto fieldY = file.readField<Scalar>(1, 1, soln.id, 1, grid_->nCells(),
                                         entry.first + "Y");

    if (fieldX.data.size() == entry.second->size() &&
        fieldY.data.size() == entry.second->size()) {
      transform(fieldX.data.begin(), fieldX.data.end(), fieldY.data.begin(),
                entry.second->begin(),
                [](Scalar x, Scalar y) { return Vector2D(x, y); });
    }

    grid_->sendMessages(*entry.second);
    entry.second->interpolateFaces();
  }

  return soln.id;
}

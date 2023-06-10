#include <set>
#include <unordered_map>
#include <vector>

#include <iostream>
#include <regex>

#include <boost/filesystem.hpp>
#include <boost/geometry.hpp>
#include <cgnslib.h>

typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>
    Node;

int main(int argc, char *argv[]) {
  using namespace std;
  using namespace boost::filesystem;

  //- Build directory lists
  regex re("Proc([0-9]+)"); // Check a processor directory
  set<path, std::function<bool(const path &, const path &)>> gridDirs(
      [re](const path &p1, const path &p2) {
        smatch matches[2];
        std::regex_search(p1.string(), matches[0], re);
        std::regex_search(p2.string(), matches[1], re);
        return stoi(matches[0][1].str()) < stoi(matches[1][1].str());
      });

  //- Load the global grids first
  for (directory_iterator end, dir("./solution"); dir != end; ++dir)
    if (regex_match(dir->path().filename().string(), re))
      gridDirs.insert(dir->path());

  //- Load transient solution data
  re = regex("[0-9]+\\.[0-9]+");
  set<path, std::function<bool(const path &, const path &)>> solutionDirs(
      [re](const path &p1, const path &p2) {
        smatch matches[2];
        std::regex_search(p1.string(), matches[0], re);
        std::regex_search(p2.string(), matches[1], re);
        return stod(matches[0][0].str()) < stod(matches[1][0].str());
      });

  //- Check if stride command exists
  int stride = 1;
  for (int argno = 1; argno < argc; ++argno)
    if (strcmp(argv[argno], "-s") == 0 || strcmp(argv[argno], "--s") == 0) {
      try {
        stride = std::stoi(argv[++argno]);
      } catch (...) {
        throw std::invalid_argument("bad argument for stride.");
      }
    } else if (strcmp(argv[argno], "-h") == 0) {
      cout << "phase-reconstruct-solution\n\n"
           << "Usage: phase-reconstruct-solution [OPTION]...\n\n"
           << "\t-s|--stride\tSpecify the spacing of time steps to "
              "reconstruct\n";

      exit(0);
    }

  for (directory_iterator end, dir("./solution"); dir != end; ++dir)
    if (regex_match(dir->path().filename().string(), re))
      solutionDirs.insert(dir->path());

  // must do this since the directory_iterator won't sort
  auto tmp = solutionDirs;
  solutionDirs.clear();

  int timeStepNo = 0;
  for (const path &p : tmp)
    if ((timeStepNo++ % stride) == 0)
      solutionDirs.insert(p);

  //- Read in global node/elements
  std::vector<Node> nodes;
  std::vector<cgsize_t> elements;
  cgsize_t nElements = 0;

  for (const auto &path : gridDirs) {
    int fn;
    cg_open((path / "Grid.cgns").c_str(), CG_MODE_READ, &fn);

    char name[256];
    cgsize_t sizes[3];
    cg_zone_read(fn, 1, 1, name, sizes);

    cgsize_t rmin = 1, rmax = sizes[0];
    std::vector<double> buffer[] = {std::vector<double>(sizes[0]),
                                    std::vector<double>(sizes[0])};

    cg_coord_read(fn, 1, 1, "CoordinateX", CGNS_ENUMV(RealDouble), &rmin, &rmax,
                  buffer[0].data());
    cg_coord_read(fn, 1, 1, "CoordinateY", CGNS_ENUMV(RealDouble), &rmin, &rmax,
                  buffer[1].data());

    cgsize_t offset = nodes.size();

    std::transform(buffer[0].begin(), buffer[0].end(), buffer[1].begin(),
                   std::back_inserter(nodes),
                   [](double x, double y) -> Node { return Node(x, y); });

    cgsize_t elementDataSize, parentData;
    cg_ElementDataSize(fn, 1, 1, 1, &elementDataSize);

    std::vector<cgsize_t> elementBuffer(elementDataSize);
    cg_elements_read(fn, 1, 1, 1, elementBuffer.data(), &parentData);

    for (int i = 0; i < elementBuffer.size();) {
      elements.push_back(elementBuffer[i]);
      switch (elementBuffer[i++]) {
      case CGNS_ENUMV(TRI_3):
        for (int j = 0; j < 3; ++j)
          elements.push_back(elementBuffer[i++] + offset);
        break;
      case CGNS_ENUMV(QUAD_4):
        for (int j = 0; j < 4; ++j)
          elements.push_back(elementBuffer[i++] + offset);
        break;
      }

      ++nElements;
    }
    cg_close(fn);
  }

  //- Compute merges
  std::vector<cgsize_t> nodeIds(nodes.size(), -1);

  {
    boost::geometry::index::rtree<pair<Node, cgsize_t>,
                                  boost::geometry::index::linear<16, 4>>
        rtree;

    cgsize_t id = 1;
    for (const Node &node : nodes)
      rtree.insert(std::make_pair(node, id++));

    double tolerance = 1e-15;
    std::vector<Node> tmp;

    id = 1;
    for (const auto &node : rtree) {
      if (nodeIds[node.second - 1] != -1) // Node already has id
        continue;

      boost::geometry::model::box<Node> box(
          Node(node.first.get<0>() - tolerance,
               node.first.get<1>() - tolerance),
          Node(node.first.get<0>() + tolerance,
               node.first.get<1>() + tolerance));

      std::vector<std::pair<Node, cgsize_t>> mergeToNode;
      rtree.query(boost::geometry::index::covered_by(box),
                  std::back_inserter(mergeToNode));

      for (const auto &mNode : mergeToNode)
        nodeIds[mNode.second - 1] = id;

      tmp.push_back(node.first);
      ++id;
    }

    std::cout << "Number of nodes: " << nodes.size() << std::endl
              << "Number of unique nodes: " << tmp.size() << std::endl;

    nodes = tmp;
  }

  std::vector<cgsize_t> elementIds(nElements, -1);

  {
    std::map<std::set<cgsize_t>, cgsize_t> elementSet;
    std::vector<cgsize_t> tmp;
    cgsize_t id = 1;
    cgsize_t nUniqueElements = 0;
    for (int i = 0; i < elements.size();) {
      std::vector<cgsize_t> element(1, elements[i]);

      switch (elements[i++]) {
      case CGNS_ENUMV(TRI_3):
        for (int j = 0; j < 3; ++j) {
          elements[i] = nodeIds[elements[i] - 1];
          element.push_back(elements[i++]);
        }
        break;
      case CGNS_ENUMV(QUAD_4):
        for (int j = 0; j < 4; ++j) {
          elements[i] = nodeIds[elements[i] - 1];
          element.push_back(elements[i++]);
        }
        break;
      }

      auto insert = elementSet.insert(std::make_pair(
          std::set<cgsize_t>(element.begin() + 1, element.end()), id));

      if (insert.second) {
        elementIds[id++ - 1] = ++nUniqueElements;
        tmp.insert(tmp.end(), element.begin(), element.end());
      } else {
        elementIds[id++ - 1] = elementIds[insert.first->second - 1];
      }
    }

    std::cout << "Number of elements: " << nElements << std::endl
              << "Number of unique elements: " << nUniqueElements << std::endl;

    nElements = nUniqueElements;
    elements = tmp;
  }

  //- Write the grid
  int fn, bid, zid, xid, sid;
  cg_open("solution.cgns", CG_MODE_WRITE, &fn);
  cg_base_write(fn, "Base", 2, 2, &bid);

  cgsize_t sizes[] = {(cgsize_t)nodes.size(), nElements, 0};
  cg_zone_write(fn, bid, "Cells", sizes, CGNS_ENUMT(Unstructured), &zid);

  std::vector<double> buffer(nodes.size());

  std::transform(nodes.begin(), nodes.end(), buffer.begin(),
                 [](const Node &node) { return node.get<0>(); });
  cg_coord_write(fn, bid, zid, CGNS_ENUMV(RealDouble), "CoordinateX",
                 buffer.data(), &xid);

  std::transform(nodes.begin(), nodes.end(), buffer.begin(),
                 [](const Node &node) { return node.get<1>(); });
  cg_coord_write(fn, bid, zid, CGNS_ENUMV(RealDouble), "CoordinateY",
                 buffer.data(), &xid);

  cg_section_write(fn, bid, zid, "Cells", CGNS_ENUMV(MIXED), 1, sizes[1],
                   sizes[2], elements.data(), &sid);

  //- Write solutions
  std::vector<double> timeValues;
  std::ostringstream flowSolutionPtrs;

  int solutionNo = 1;

  for (const auto &path : solutionDirs) {
    double time = std::stod(path.filename().string());
    timeValues.push_back(time);

    std::unordered_map<std::string, std::vector<double>> doubleFields;
    std::unordered_map<std::string, std::vector<int>> intFields;

    for (const auto &subpath : gridDirs) {
      std::cout << "Loading: " << path / subpath.filename() / "Solution.cgns"
                << "..." << std::endl;

      int fn2;
      cg_open((path / subpath.filename() / "Solution.cgns").c_str(),
              CG_MODE_READ, &fn2);

      char name[256];
      cg_zone_read(fn2, 1, 1, name, sizes);

      int nfields = 0;
      cg_nfields(fn2, 1, 1, 1, &nfields);

      for (int field = 1; field <= nfields; ++field) {
        CGNS_ENUMT(DataType_t) type;
        cg_field_info(fn2, 1, 1, 1, field, &type, name);
        cgsize_t rmin = 1, rmax = sizes[1];

        switch (type) {
        case CGNS_ENUMV(Integer): {
          std::vector<int> buffer(sizes[1]);
          cg_field_read(fn2, 1, 1, 1, name, CGNS_ENUMV(Integer), &rmin, &rmax,
                        buffer.data());
          intFields[name].insert(intFields[name].end(), buffer.begin(),
                                 buffer.end());
        } break;
        case CGNS_ENUMV(RealDouble): {
          std::vector<double> buffer(sizes[1]);
          cg_field_read(fn2, 1, 1, 1, name, CGNS_ENUMV(RealDouble), &rmin,
                        &rmax, buffer.data());
          doubleFields[name].insert(doubleFields[name].end(), buffer.begin(),
                                    buffer.end());
        } break;
        }
      }

      cg_close(fn2);
    }

    std::string flowSolutionPtr = "FlowSolution" + std::to_string(solutionNo);
    flowSolutionPtrs << setw(32) << setfill(' ') << left << flowSolutionPtr;

    int sid;
    cg_sol_write(fn, 1, 1, flowSolutionPtr.c_str(), CGNS_ENUMV(CellCenter),
                 &sid);

    for (const auto &field : intFields) {
      int fieldId;
      std::vector<int> buffer(nElements);
      for (int i = 0; i < field.second.size(); ++i)
        buffer[elementIds[i] - 1] = field.second[i];

      cg_field_write(fn, 1, 1, sid, CGNS_ENUMV(Integer), field.first.c_str(),
                     buffer.data(), &fieldId);
    }

    for (const auto &field : doubleFields) {
      int fieldId;
      std::vector<double> buffer(nElements);
      for (int i = 0; i < field.second.size(); ++i)
        buffer[elementIds[i] - 1] = field.second[i];

      cg_field_write(fn, 1, 1, sid, CGNS_ENUMV(RealDouble), field.first.c_str(),
                     buffer.data(), &fieldId);
    }

    ++solutionNo;
  }

  //- Write zone iterative data
  cg_ziter_write(fn, 1, 1, "ZoneIterativeData");
  cg_goto(fn, 1, "Zone_t", 1, "ZoneIterativeData_t", 1, "end");
  sizes[0] = 32;
  sizes[1] = timeValues.size();
  cg_array_write("FlowSolutionPointers", CGNS_ENUMV(Character), 2, sizes,
                 flowSolutionPtrs.str().c_str());

  //- Write base iterative data
  cg_biter_write(fn, 1, "TimeIterValues", timeValues.size());
  cg_goto(fn, 1, "BaseIterativeData_t", 1, "end");
  cg_array_write("TimeValues", CGNS_ENUMV(RealDouble), 1, &sizes[1],
                 timeValues.data());
  cg_simulation_type_write(fn, 1, CGNS_ENUMV(TimeAccurate));

  //- Finalize
  cg_close(fn);

  return 0;
}

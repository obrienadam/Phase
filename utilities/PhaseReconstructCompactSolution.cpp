#include <iostream>
#include <numeric>
#include <regex>
#include <set>
#include <unordered_map>

#include <boost/filesystem.hpp>
#include <boost/geometry.hpp>
#include <cgnslib.h>

typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>
    Node;
typedef boost::geometry::index::rtree<std::pair<Node, size_t>,
                                      boost::geometry::index::quadratic<16, 4>>
    RTree;

int main(int argc, char *argv[]) {
  using namespace std;
  using namespace boost::filesystem;
  namespace bgi = boost::geometry::index;

  //- Collect a list of all the output files
  regex re("proc(\\d+)\\.cgns");
  set<path, function<bool(const path &, const path &)>> files(
      [re](const path &p1, const path &p2) {
        smatch matches[2];
        regex_search(p1.string(), matches[0], re);
        regex_search(p2.string(), matches[1], re);
        return stoi(matches[0][1].str()) < stoi(matches[1][1].str());
      });

  for (directory_iterator end, dir("./solution"); dir != end; ++dir)
    if (regex_match(dir->path().filename().string(), re))
      files.insert(dir->path());

  //- global nodes and elements
  vector<double> buffer, xcoords, ycoords;
  vector<int> cptr(1, 0), cind;
  vector<cgsize_t> elements;

  //- ownerhship and globalIds
  vector<int> ownership, globalIds;

  //- Flow solution points
  vector<array<char, 32>> flowSolutionPtrs;
  vector<double> timeValues;

  //- File ids, so files don't need to be open and closed repeatedly (this seems
  //to cause major performance issues)
  vector<int> fids;

  //- Read the files
  int proc = 0, nodeIdStart = 0;
  for (const path &p : files) {
    //- Do a check to make sure file corresponds to the expected proc
    smatch match;
    regex_search(p.string(), match, re);

    if (stoi(match[1]) != proc)
      throw runtime_error("Missing grid file for proc " + to_string(proc) +
                          ".");

    cout << "Reading file " << p.c_str() << "...\n";

    int fid;
    cg_open(p.c_str(), CG_MODE_READ, &fid);
    fids.emplace_back(fid);

    array<char, 256> name;
    int cellDim, physDim;

    cg_base_read(fid, 1, name.data(), &cellDim, &physDim);

    if (cellDim != 2 || physDim != 2)
      throw runtime_error(string("Bad cell or physical dimension for base \"") +
                          name.data() + "\".");

    cout << "\t--Read base \"" << name.data() << "\".\n";

    cgsize_t sizes[3];
    cg_zone_read(fid, 1, 1, name.data(), sizes);

    cout << "\t--Read zone \"" << name.data() << "\".\n"
         << "\t\t--Num nodes = " << sizes[0] << "\n"
         << "\t\t--Num cells = " << sizes[1] << "\n";

    cgsize_t rmin = 1, rmax = sizes[0];
    buffer.resize(sizes[0]);

    cg_coord_read(fid, 1, 1, "CoordinateX", CGNS_ENUMV(RealDouble), &rmin,
                  &rmax, buffer.data());
    xcoords.insert(xcoords.end(), buffer.begin(), buffer.end());

    cg_coord_read(fid, 1, 1, "CoordinateY", CGNS_ENUMV(RealDouble), &rmin,
                  &rmax, buffer.data());
    ycoords.insert(ycoords.end(), buffer.begin(), buffer.end());

    CGNS_ENUMT(GridLocation_t) location;

    cg_sol_info(fid, 1, 1, 1, name.data(), &location);

    if (name.data() != string("Info")) {
      cerr << "First solution node must contain fields \"ProcNo\" and "
              "\"GlobalID\" and be named \"Info\".\n";
      throw exception();
    }

    //- Get the element ownership data and ids
    ownership.resize(sizes[1]);
    globalIds.resize(sizes[1]);

    rmin = 1;
    rmax = sizes[1];
    cg_field_read(fid, 1, 1, 1, "ProcNo", CGNS_ENUMV(Integer), &rmin, &rmax,
                  ownership.data());
    cg_field_read(fid, 1, 1, 1, "GlobalID", CGNS_ENUMV(Integer), &rmin, &rmax,
                  globalIds.data());

    cgsize_t elemDataSize;
    cg_ElementDataSize(fid, 1, 1, 1, &elemDataSize);

    elements.resize(elemDataSize);
    cg_elements_read(fid, 1, 1, 1, elements.data(), nullptr);

    for (size_t i = 0, elemNo = 0; i < elements.size(); ++elemNo) {
      size_t npts;
      switch (elements[i++]) {
      case CGNS_ENUMV(TRI_3):
        npts = 3;
        break;
      case CGNS_ENUMV(QUAD_4):
        npts = 4;
        break;
      default:
        cerr << "Only TRI_3 and QUAD_4 element types are supported.\n";
        throw exception();
      }

      //- Don't add elements that are not owned by this proc
      if (ownership[elemNo] == proc) {
        for_each(elements.begin() + i, elements.begin() + i + npts,
                 [nodeIdStart](cgsize_t &id) { id += nodeIdStart - 1; });
        cptr.emplace_back(cptr.back() + npts);
        cind.insert(cind.end(), elements.begin() + i,
                    elements.begin() + i + npts);
      }

      i += npts;
    }

    //- Construct the flow solution pointers
    int nsols;
    cg_nsols(fid, 1, 1, &nsols);

    bool populateFlowSolutionPointers = flowSolutionPtrs.empty();

    for (int S = 2; S <= nsols; ++S) {
      array<char, 32> flowSolutionPtr = {'\0'};
      CGNS_ENUMT(GridLocation_t) location;
      cg_sol_info(fid, 1, 1, S, flowSolutionPtr.data(), &location);

      //- Make sure is is the expected flow solution pointer
      if ("FlowSolution" + to_string(S - 1) != flowSolutionPtr.data())
        throw runtime_error("Expected flow solution pointer \"FlowSolution" +
                            to_string(S - 1) + "\", but received \"" +
                            flowSolutionPtr.data() + "\".");

      //- First proc is used to deduce the number of flow solution pointers as
      //well as the time values

      if (populateFlowSolutionPointers) {
        flowSolutionPtrs.emplace_back(flowSolutionPtr);
        char *text;
        cg_goto(fid, 1, "Zone_t", 1, "FlowSolution_t", S, "end");
        cg_descriptor_read(1, name.data(), &text);
        timeValues.emplace_back(stod(text));
        cg_free(text);
      }
    }

    if (nsols - 2 + 1 < flowSolutionPtrs.size())
      throw runtime_error("Too few flow solution pointers on proc " +
                          to_string(proc) + ".");
    else if (nsols - 2 + 1 > flowSolutionPtrs.size())
      cout << "\t*** WARNING *** Additional flow solution pointers on proc "
           << proc << ". Expected " << flowSolutionPtrs.size() << ", received "
           << nsols - 2 + 1 << ".\n";

    proc++;
    nodeIdStart += sizes[0];
  }

  //- Perform node merging between procs

  cout << "Construction of global cells complete. Number of global cells = "
       << cptr.size() - 1 << ".\n"
       << "Merging duplicate nodes...\n";

  RTree rtree;
  for (size_t i = 0; i < xcoords.size(); ++i)
    rtree.insert(make_pair(Node(xcoords[i], ycoords[i]), i));

  const double tolerance = 1e-13;
  vector<size_t> merges(rtree.size());
  iota(merges.begin(), merges.end(), 0);

  size_t nMerges = 0;
  for (const RTree::value_type &v : rtree) {
    //- Check if node is already being merged, no need to check again if found
    if (merges[v.second] != v.second)
      continue;

    //- Construct search box around node
    boost::geometry::model::box<Node> box(
        Node(v.first.get<0>() - tolerance, v.first.get<1>() - tolerance),
        Node(v.first.get<0>() + tolerance, v.first.get<1>() + tolerance));

    auto p = bgi::covered_by(
        box); //- merge all nodes at this location into this node
    for (auto qit = rtree.qbegin(p); qit != rtree.qend(); ++qit) {
      merges[qit->second] = v.second;
      if (qit->second != v.second)
        ++nMerges;
    }
  }
  cout << "Merging " << nMerges << " duplicate nodes.\n";

  vector<double> newXCoords, newYCoords;
  vector<size_t> newNodeIds(merges.size());
  for (size_t i = 0, id = 0; i < merges.size(); ++i)
    if (merges[i] == i) {
      newXCoords.emplace_back(xcoords[i]);
      newYCoords.emplace_back(ycoords[i]);
      newNodeIds[i] = id++;
    }

  for (size_t i = 0; i < merges.size(); ++i)
    if (merges[i] != i)
      newNodeIds[i] = newNodeIds[merges[i]];

  xcoords = move(newXCoords);
  ycoords = move(newYCoords);

  for (size_t i = 0; i < cind.size(); ++i)
    cind[i] = newNodeIds[cind[i]];

  //- Create new elements data structure for writing
  elements.clear();
  for (size_t i = 0; i < cptr.size() - 1; ++i) {
    switch (cptr[i + 1] - cptr[i]) {
    case 3:
      elements.emplace_back(CGNS_ENUMV(TRI_3));
      break;
    case 4:
      elements.emplace_back(CGNS_ENUMV(QUAD_4));
      break;
    }

    for (int j = cptr[i]; j < cptr[i + 1]; ++j)
      elements.emplace_back(cind[j] + 1);
  }

  //- Write the grid info to the output file

  int fid;
  cg_open("solution.cgns", CG_MODE_WRITE, &fid);

  int bid;
  cg_base_write(fid, "solution", 2, 2, &bid);
  cg_simulation_type_write(fid, 1, CGNS_ENUMV(TimeAccurate));

  int zid;
  cgsize_t sizes[3] = {cgsize_t(xcoords.size()), cgsize_t(cptr.size() - 1), 0};
  cg_zone_write(fid, bid, "Cells", sizes, CGNS_ENUMT(Unstructured), &zid);

  int cid;
  cg_coord_write(fid, bid, zid, CGNS_ENUMV(RealDouble), "CoordinateX",
                 xcoords.data(), &cid);
  cg_coord_write(fid, bid, zid, CGNS_ENUMV(RealDouble), "CoordinateY",
                 ycoords.data(), &cid);

  int sid;
  cg_section_write(fid, bid, zid, "Elements", CGNS_ENUMV(MIXED), 1,
                   cptr.size() - 1, 0, elements.data(), &sid);

  //- write all fields
  vector<int> i_buffers[2];
  vector<double> d_buffers[2];

  //- Create transient flow solution nodes
  for (const auto &flowSolutionPtr : flowSolutionPtrs)
    cg_sol_write(fid, 1, 1, flowSolutionPtr.data(), CGNS_ENUMV(CellCenter),
                 &sid);

  //- Create grid info solution node
  cg_sol_write(fid, 1, 1, "Info", CGNS_ENUMV(CellCenter), &sid);

  cgsize_t rglobal[2] = {1, 1};

  for (int proc = 0; proc < fids.size(); ++proc) {
    cout << "Writing data for proc " << proc << "...\n";

    int fid2 = fids[proc];

    char name[256];
    cgsize_t sizes[3];
    cg_zone_read(fid2, 1, 1, name, sizes);

    cgsize_t rlocal[2] = {1, sizes[1]};

    ownership.resize(sizes[1]);
    cg_field_read(fid2, 1, 1, 1, "ProcNo", CGNS_ENUMV(Integer), &rlocal[0],
                  &rlocal[1], ownership.data());

    int ncells = count_if(ownership.begin(), ownership.end(),
                          [proc](int p) { return p == proc; });

    //- Get the upper bound for all writes
    rglobal[1] = rglobal[0] + ncells - 1;

    //- Write ownership data
    i_buffers[0].resize(rlocal[1]);
    i_buffers[1].resize(rglobal[1] - rglobal[0] + 1);

    i_buffers[0] = ownership;
    for (int i = 0, j = 0; i < ownership.size(); ++i)
      if (ownership[i] == proc)
        i_buffers[1][j++] = i_buffers[0][i];

    int F;
    cg_field_partial_write(fid, 1, 1, flowSolutionPtrs.size() + 1,
                           CGNS_ENUMV(Integer), "ProcNo", &rglobal[0],
                           &rglobal[1], i_buffers[1].data(), &F);

    cg_field_read(fid2, 1, 1, 1, "GlobalID", CGNS_ENUMV(Integer), &rlocal[0],
                  &rlocal[1], i_buffers[0].data());

    for (int i = 0, j = 0; i < ownership.size(); ++i)
      if (ownership[i] == proc)
        i_buffers[1][j++] = i_buffers[0][i];

    cg_field_partial_write(fid, 1, 1, flowSolutionPtrs.size() + 1,
                           CGNS_ENUMV(Integer), "GlobalID", &rglobal[0],
                           &rglobal[1], i_buffers[1].data(), &F);

    //- Write flow solution data
    d_buffers[0].resize(rlocal[1]);
    d_buffers[1].resize(rglobal[1] - rglobal[0] + 1);

    for (int S = 1; S <= flowSolutionPtrs.size(); ++S) {
      int nfields;
      cg_nfields(fid2, 1, 1, S + 1, &nfields);

      for (int F = 1; F <= nfields; ++F) {
        array<char, 256> name;
        CGNS_ENUMT(DataType_t) type;
        cg_field_info(fid2, 1, 1, S + 1, F, &type, name.data());
        cg_field_read(fid2, 1, 1, S + 1, name.data(), CGNS_ENUMV(RealDouble),
                      &rlocal[0], &rlocal[1], d_buffers[0].data());

        for (int i = 0, j = 0; i < ownership.size(); ++i)
          if (ownership[i] == proc)
            d_buffers[1][j++] = d_buffers[0][i];

        int F2;
        cg_field_partial_write(fid, 1, 1, S, CGNS_ENUMV(RealDouble),
                               name.data(), &rglobal[0], &rglobal[1],
                               d_buffers[1].data(), &F2);
      }
    }

    //- File is no longer needed and can be closed
    cg_close(fid2);
    rglobal[0] = rglobal[1] + 1;
  }

  //- Write base iterative data
  cg_biter_write(fid, 1, "TimeIterValues", flowSolutionPtrs.size());
  cg_goto(fid, 1, "BaseIterativeData_t", 1, "end");

  sizes[0] = (cgsize_t)timeValues.size();
  cg_array_write("TimeValues", CGNS_ENUMV(RealDouble), 1, sizes,
                 timeValues.data());

  //- Write sone iterative data
  cg_ziter_write(fid, 1, 1, "ZoneIterativeData");
  cg_goto(fid, 1, "Zone_t", 1, "ZoneIterativeData_t", 1, "end");

  sizes[0] = 32;
  sizes[1] = flowSolutionPtrs.size();

  cg_array_write("FlowSolutionPointers", CGNS_ENUMV(Character), 2, sizes,
                 flowSolutionPtrs.data());
  cg_close(fid);
}

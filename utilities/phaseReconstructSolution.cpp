#include <vector>
#include <regex>
#include <iostream>

#include <cgnslib.h>
#include <boost/geometry.hpp>
#include <boost/filesystem.hpp>

typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> Point;

struct Section
{
    char name[32];
    std::vector<cgsize_t> elements;
    CGNS_ENUMT(ElementType_t) type;
    cgsize_t start, end, parentData;
    int nbndry, parentFlag;
};

struct Field
{
    char name[32];
    std::vector<double> scalarData;
    std::vector<int> integerData;
    CGNS_ENUMT(DataType_t) type;
};

struct Solution
{
    char name[32];
    std::vector<Field> fields;
    CGNS_ENUMT(GridLocation_t) location;
};

struct Zone
{
    char name[32];
    cgsize_t sizes[3];
    std::vector<Point> nodes;
    std::vector<Section> sections;
    Solution solution;
    CGNS_ENUMT(ZoneType_t) type;
};

struct Base
{
    char name[32];
    int cellDim, physDim;
    std::vector<Zone> zones;
};

Base loadGrid(const std::string &filename);

Solution loadSolution(const std::string &filename);

void mergeGrids(const std::vector<Base> &grids,
                const std::vector<std::vector<Solution>> &solutions,
                std::vector<Point> &nodes, std::vector<cgsize_t> &elements, int nTimeSteps,
                std::map<std::string, std::vector<Field>> &fields);

std::vector<cgsize_t> getNodeMergeList(const std::vector<Point> &nodes);

std::vector<Point> mergeNodes(const std::vector<Point> &nodes,
                              const std::vector<cgsize_t> &nodeMergeList,
                              std::vector<cgsize_t> &elements);

std::vector<cgsize_t> mergeElements(cgsize_t nNodes,
                                    const std::vector<cgsize_t> &elements,
                                    int &nElements,
                                    std::map<std::string, std::vector<Field> > &fields);

int main(int argc, char *argv[])
{
    using namespace std;
    using namespace boost::filesystem;

    regex re("Proc[0-9]+"); // Check a processor directory
    set<path> gridDirs;

    //- Load the global grids first
    for (directory_iterator end, dir("./solution"); dir != end; ++dir)
        if (regex_match(dir->path().filename().c_str(), re))
        {
            path gridFile = dir->path();
            gridDirs.insert(gridFile /= "Grid.cgns");
        }

    vector<Base> procGrids;
    for (path p: gridDirs)
        procGrids.push_back(loadGrid(p.c_str()));

    //- Load transient solution data
    re = regex("[0-9]+\\.[0-9]+");

    std::vector<path> solutionDirs;
    for (directory_iterator end, dir("./solution"); dir != end; ++dir)
        if (regex_match(dir->path().filename().c_str(), re))
            solutionDirs.push_back(dir->path());

    std::sort(solutionDirs.begin(), solutionDirs.end(), [](const path& p1, const path& p2){
        double t1 = stod(p1.filename().c_str());
        double t2 = stod(p2.filename().c_str());

        return t1 < t2;
    });

    vector<double> timeSteps;
    vector<vector<Solution>> procSolutions(procGrids.size());

    for (path solutionDir: solutionDirs)
    {
        int proc = 0;
        for (path procDir: gridDirs)
        {
            procDir = solutionDir / procDir.parent_path().filename() / "Solution.cgns";
            procSolutions[proc++].push_back(loadSolution(procDir.c_str()));
        }

        timeSteps.push_back(stod(solutionDir.filename().c_str()));
    }

    vector<Point> nodes;
    vector<cgsize_t> elements;
    map<string, vector<Field>> fields;

    mergeGrids(procGrids, procSolutions, nodes, elements, timeSteps.size(), fields);
    nodes = mergeNodes(nodes, getNodeMergeList(nodes), elements);

    int nElements;
    elements = mergeElements(nodes.size(), elements, nElements, fields);

    cout << "Writing a new .cgns file...\n";

    //- Write a new merged cgns file

    cgsize_t sizes[3];
    sizes[0] = nodes.size();
    sizes[1] = nElements;
    sizes[2] = 0;

    int fid;
    cg_open("solution.cgns", CG_MODE_WRITE, &fid);

    int bid;
    cg_base_write(fid, "Base", 2, 2, &bid);

    int zid;
    cg_zone_write(fid, bid, "Cells", sizes, CGNS_ENUMT(Unstructured), &zid);

    vector<double> xCoords(sizes[0]), yCoords(sizes[0]);
    std::transform(nodes.begin(), nodes.end(), xCoords.begin(), [](const Point &pt) { return pt.get<0>(); });
    std::transform(nodes.begin(), nodes.end(), yCoords.begin(), [](const Point &pt) { return pt.get<1>(); });

    int xid;
    cg_coord_write(fid, bid, zid, CGNS_ENUMV(RealDouble), "CoordinateX", xCoords.data(), &xid);
    cg_coord_write(fid, bid, zid, CGNS_ENUMV(RealDouble), "CoordinateY", yCoords.data(), &xid);

    int sid;
    cg_section_write(fid, bid, zid, "Cells", CGNS_ENUMV(MIXED), 1, sizes[1], 0, elements.data(), &sid);

    ostringstream sout;

    for (int i = 0; i < timeSteps.size(); ++i)
    {
        int solId;

        string flowSolutionPointer = "FlowSolution" + to_string(i + 1);

        cg_sol_write(fid, bid, zid, flowSolutionPointer.c_str(), CGNS_ENUMV(CellCenter), &solId);

        for (const auto &field: fields)
        {
            int fieldId;

            switch (field.second[i].type)
            {
                case CGNS_ENUMV(Integer):
                    cg_field_write(fid, bid, zid, solId, CGNS_ENUMV(Integer), field.first.c_str(),
                                   field.second[i].integerData.data(), &fieldId);
                    break;
                case CGNS_ENUMV(RealDouble):
                    cg_field_write(fid, bid, zid, solId, CGNS_ENUMV(RealDouble), field.first.c_str(),
                                   field.second[i].scalarData.data(), &fieldId);
                    break;
            }
        }

        sout.width(32);
        sout << left << setfill(' ') << flowSolutionPointer;
    }

    //- Write zone iterative data
    cg_ziter_write(fid, bid, zid, "ZoneIterativeData");
    cg_goto(fid, bid, "Zone_t", zid, "ZoneIterativeData_t", 1, "end");

    sizes[0] = 32;
    sizes[1] = timeSteps.size();

    cg_array_write("FlowSolutionPointers", CGNS_ENUMV(Character), 2, sizes, sout.str().c_str());

    //- Write base iterative data
    cg_biter_write(fid, bid, "TimeIterValues", timeSteps.size());
    cg_goto(fid, bid, "BaseIterativeData_t", 1, "end");

    cg_array_write("TimeValues", CGNS_ENUMV(RealDouble), 1, &sizes[1], timeSteps.data());

    cg_simulation_type_write(fid, bid, CGNS_ENUMV(TimeAccurate));

    cg_close(fid);

    cout << "Solution successfully reconstructed.\n";

    return 0;
}

Base loadGrid(const std::string &filename)
{
    int fid;

    std::cout << "Loading grid \"" << filename << "\"...\n";

    cg_open(filename.c_str(), CG_MODE_READ, &fid);

    Base base;

    cg_base_read(fid, 1, base.name, &base.cellDim, &base.physDim);

    int nZones;
    cg_nzones(fid, 1, &nZones);

    for (int zid = 1; zid <= nZones; ++zid)
    {
        base.zones.push_back(Zone());
        Zone &zone = base.zones.back();

        cg_zone_read(fid, 1, zid, zone.name, zone.sizes);
        cg_zone_type(fid, 1, zid, &zone.type);

        std::vector<double> xCoords(zone.sizes[0]), yCoords(zone.sizes[0]);

        cgsize_t rmin = 1, rmax = zone.sizes[0];
        cg_coord_read(fid, 1, zid, "CoordinateX", CGNS_ENUMV(RealDouble), &rmin, &rmax, xCoords.data());
        cg_coord_read(fid, 1, zid, "CoordinateY", CGNS_ENUMV(RealDouble), &rmin, &rmax, yCoords.data());

        zone.nodes.resize(zone.sizes[0]);
        std::transform(xCoords.begin(), xCoords.end(), yCoords.begin(), zone.nodes.begin(),
                       [](double x, double y) { return Point(x, y); });

        int nSections;
        cg_nsections(fid, 1, zid, &nSections);

        for (int sid = 1; sid <= nSections; ++sid)
        {
            zone.sections.push_back(Section());
            Section &section = zone.sections.back();

            cg_section_read(fid, 1, zid, sid, section.name, &section.type, &section.start, &section.end,
                            &section.nbndry, &section.parentFlag);

            cgsize_t elementDataSize;
            cg_ElementDataSize(fid, 1, zid, sid, &elementDataSize);

            section.elements.resize(elementDataSize);
            cg_elements_read(fid, 1, zid, sid, section.elements.data(), &section.parentData);
        }

        int nbocos;
        cg_nbocos(fid, 1, zid, &nbocos);
    }

    cg_close(fid);

    return base;
}

Solution loadSolution(const std::string &filename)
{
    Solution solution;
    int fid;

    std::cout << "Loading solution \"" << filename << "\"...\n";

    cg_open(filename.c_str(), CG_MODE_READ, &fid);

    cg_sol_info(fid, 1, 1, 1, solution.name, &solution.location);

    int datadim;
    cgsize_t dimvals[3];
    cg_sol_size(fid, 1, 1, 1, &datadim, dimvals);

    if (datadim != 1)
    {
        std::cout << "Solution data dimension is not 1, exiting...\n";
        exit(-1);
    }

    cgsize_t rmin = 1, rmax = dimvals[0];

    int nfields;
    cg_nfields(fid, 1, 1, 1, &nfields);

    for (int fieldid = 1; fieldid <= nfields; ++fieldid)
    {
        solution.fields.push_back(Field());
        Field &field = solution.fields.back();

        cg_field_info(fid, 1, 1, 1, fieldid, &field.type, field.name);

        switch (field.type)
        {
            case CGNS_ENUMV(Integer):
                field.integerData.resize(dimvals[0]);
                cg_field_read(fid, 1, 1, 1, field.name, field.type, &rmin, &rmax, field.integerData.data());
                break;
            case CGNS_ENUMV(RealDouble):
                field.scalarData.resize(dimvals[0]);
                cg_field_read(fid, 1, 1, 1, field.name, field.type, &rmin, &rmax, field.scalarData.data());
                break;
            default:
                std::cout << "Urecognized datat type.\n";
                exit(-1);
        }
    }

    cg_close(fid);

    return solution;
}

void mergeGrids(const std::vector<Base> &grids,
                const std::vector<std::vector<Solution>> &solutions,
                std::vector<Point> &nodes, std::vector<cgsize_t> &elements,
                int nTimeSteps,
                std::map<std::string, std::vector<Field> > &fields)
{
    using namespace std;

    nodes.clear();
    elements.clear();
    fields.clear();

    int gridNo = 0;
    for (const Base &grid: grids)
        for (const Zone &zone: grid.zones)
        {
            cgsize_t offset = nodes.size();
            nodes.insert(nodes.end(), zone.nodes.begin(), zone.nodes.end());

            for (const Section &section: zone.sections)
            {
                if (section.type != CGNS_ENUMV(MIXED))
                    continue;

                for (cgsize_t i = 0; i < section.elements.size();)
                {
                    cgsize_t type = section.elements[i++];
                    elements.push_back(type);

                    switch (type)
                    {
                        case CGNS_ENUMV(TRI_3):

                            for (cgsize_t j = 0; j < 3; ++j, ++i)
                                elements.push_back(section.elements[i] + offset);

                            break;
                        case CGNS_ENUMV(QUAD_4):
                            for (cgsize_t j = 0; j < 4; ++j, ++i)
                                elements.push_back(section.elements[i] + offset);

                            break;
                        default:
                        {
                            cout << "Invalid element type. Exiting...\n";
                            exit(-1);
                        }
                    };
                }

                int timeStepNo = 0;
                for (const Solution &soln: solutions[gridNo++])
                {
                    for (const Field &field: soln.fields)
                    {
                        vector<Field> &globalField = fields[field.name];

                        if (globalField.size() < nTimeSteps)
                            globalField.resize(nTimeSteps);

                        globalField[timeStepNo].type = field.type;

                        switch (field.type)
                        {
                            case CGNS_ENUMV(Integer):
                                globalField[timeStepNo].integerData.insert(globalField[timeStepNo].integerData.end(),
                                                                           field.integerData.begin(),
                                                                           field.integerData.end());
                                break;
                            case CGNS_ENUMV(RealDouble):
                                globalField[timeStepNo].scalarData.insert(globalField[timeStepNo].scalarData.end(),
                                                                          field.scalarData.begin(),
                                                                          field.scalarData.end());
                                break;
                        }
                    }

                    ++timeStepNo;
                }
            }
        }
}

std::vector<cgsize_t> getNodeMergeList(const std::vector<Point> &nodes)
{
    using namespace std;

    double tolerance = 1e-15;

    cout << "Searching " << nodes.size() << " nodes for duplicates, tolerance = " << tolerance << "...\n";

    //- Construct a searchable tree
    boost::geometry::index::rtree<pair<Point, cgsize_t>, boost::geometry::index::quadratic<32>> rtree;

    cgsize_t id = 1;
    for (const Point &node: nodes)
        rtree.insert(make_pair(node, id++));

    //- Find nodes that are the same within some tolerance, mark them for merging
    vector<cgsize_t> nodeMerges(nodes.size(), 0);
    int nDuplicates = 0;

    id = 0;
    for (const Point &node: nodes)
    {
        ++id;

        if (nodeMerges[id - 1] != 0) // Node is already marked for merging
            continue;

        boost::geometry::model::box<Point> box(
                Point(node.get<0>() - tolerance, node.get<1>() - tolerance),
                Point(node.get<0>() + tolerance, node.get<1>() + tolerance)
        );

        std::vector<pair<Point, cgsize_t>> result;
        rtree.query(boost::geometry::index::covered_by(box), back_inserter(result));

        for (const auto &val: result)
        {
            if (val.second == id)
                continue;

            if (nodeMerges[val.second - 1] != 0)
                continue;

            nodeMerges[val.second - 1] = id; // Works because the rtree returns all the nodes at this location

            ++nDuplicates;
        }
    }

    cout << "Found " << nDuplicates << " duplicate nodes.\n";

    return nodeMerges;
}

std::vector<Point>
mergeNodes(const std::vector<Point> &nodes, const std::vector<cgsize_t> &nodeMergeList, std::vector<cgsize_t> &elements)
{
    //- Merge the nodes that have been marked. Update the elements list accordingly
    std::vector<Point> mergedNodes;
    std::vector<cgsize_t> mergedNodeIds(nodes.size(), 0);

    cgsize_t id = 1;
    for (int i = 0; i < nodes.size(); ++i)
        if (nodeMergeList[i] == 0) // not merged into another node, gets a new id
        {
            mergedNodeIds[i] = id++;
            mergedNodes.push_back(nodes[i]);
        }

    for (cgsize_t i = 0; i < elements.size();)
    {
        cgsize_t type = elements[i++];

        cgsize_t nVerts;

        switch (type)
        {
            case CGNS_ENUMV(TRI_3):
                nVerts = 3;
                break;
            case CGNS_ENUMV(QUAD_4):
                nVerts = 4;
                break;
            default:
                std::cout << "Invalid element type. Exiting...\n";
                exit(-1);
        }

        for (cgsize_t j = 0; j < nVerts; ++j, ++i)
        {
            if (nodeMergeList[elements[i] - 1])
                elements[i] = mergedNodeIds[nodeMergeList[elements[i] - 1] - 1];
            else
                elements[i] = mergedNodeIds[elements[i] - 1];
        }
    }

    return mergedNodes;
}

std::vector<cgsize_t> mergeElements(cgsize_t nNodes,
                                    const std::vector<cgsize_t> &elements, int &nElements,
                                    std::map<std::string, std::vector<Field> > &fields)
{
    using namespace std;

    cout << "Searching for duplicate elements...\n";

    //- Merge duplicate elements
    //- Construct a node to element map
    vector<vector<cgsize_t>> nodeToElements(nNodes);
    vector<vector<cgsize_t>> elementToNodes;

    cgsize_t id = 1;
    for (cgsize_t i = 0; i < elements.size(); ++id)
    {
        cgsize_t type = elements[i++];

        cgsize_t nVerts;

        switch (type)
        {
            case CGNS_ENUMV(TRI_3):
                nVerts = 3;
                break;
            case CGNS_ENUMV(QUAD_4):
                nVerts = 4;
                break;
            default:
                cout << "Invalid element type. Exiting...\n";
                exit(-1);
        }

        elementToNodes.push_back(vector<cgsize_t>());

        for (cgsize_t j = 0; j < nVerts; ++j, ++i)
        {
            nodeToElements[elements[i] - 1].push_back(id);
            elementToNodes[id - 1].push_back(elements[i]);
        }
    }

    //- Loop over each node, check for duplicate elements and mark them for removal
    vector<bool> removeElement(elementToNodes.size(), false);
    int nDuplicates = 0;
    for (const vector<cgsize_t> &elements: nodeToElements)
    {
        for (cgsize_t ele1: elements)
        {
            if (removeElement[ele1 - 1]) // don't check elements marked for removal
                continue;

            for (cgsize_t ele2: elements)
            {
                if (ele1 == ele2 || removeElement[ele2 - 1])
                    continue;

                bool duplicate = true;

                for (cgsize_t node1: elementToNodes[ele1 - 1])
                    duplicate &= find(elementToNodes[ele2 - 1].begin(), elementToNodes[ele2 - 1].end(), node1) !=
                                 elementToNodes[ele2 - 1].end();

                if (duplicate)
                {
                    removeElement[ele2 - 1] = true;
                    ++nDuplicates;
                }
            }
        }
    }

    cout << "Found " << nDuplicates << " duplicate elements.\n";

    //- Construct new element list
    vector<cgsize_t> newElements;
    id = 1;
    nElements = 0;
    for (cgsize_t i = 0; i < elements.size(); ++id)
    {
        cgsize_t type = elements[i++];

        cgsize_t nVerts;

        switch (type)
        {
            case CGNS_ENUMV(TRI_3):
                nVerts = 3;
                break;
            case CGNS_ENUMV(QUAD_4):
                nVerts = 4;
                break;
            default:
                cout << "Invalid element type. Exiting...\n";
                exit(-1);
        }

        if (removeElement[id - 1])
        {
            i += nVerts;
            continue;
        }

        newElements.push_back(type);

        for (cgsize_t j = 0; j < nVerts; ++j, ++i)
            newElements.push_back(elements[i]);

        ++nElements;
    }

    //- Modify field data
    for (auto &field: fields)
    {
        for (Field &timeStep: field.second)
        {
            switch (timeStep.type)
            {
                case CGNS_ENUMV(Integer):
                {
                    vector<int> newTimeStep;

                    auto removeItr = removeElement.begin();

                    for (int val: timeStep.integerData)
                        if (!*(removeItr++))
                            newTimeStep.push_back(val);

                    timeStep.integerData = newTimeStep;
                }
                    break;
                case CGNS_ENUMV(RealDouble):
                {
                    vector<double> newTimeStep;

                    auto removeItr = removeElement.begin();

                    for (double val: timeStep.scalarData)
                        if (!*(removeItr++))
                            newTimeStep.push_back(val);

                    timeStep.scalarData = newTimeStep;
                }
                    break;
            }
        }
    }

    return newElements;
}

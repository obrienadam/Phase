#include "CgnsUnstructuredGrid.h"
#include "Exception.h"

CgnsUnstructuredGrid::CgnsUnstructuredGrid()
    :
      FiniteVolumeGrid2D()
{
    fileIsOpen_ = false;
}

CgnsUnstructuredGrid::CgnsUnstructuredGrid(const Input &input)
    :
      CgnsUnstructuredGrid()
{
    const std::string filename = input.caseInput().get<std::string>("Grid.filename");
    const Scalar convertToMeters = input.caseInput().get<Scalar>("Grid.convertToMeters", 1.);

    int fileId;
    char name[256];

    cg_open(filename.c_str(), CG_MODE_READ, &fileId);

    //- The first base is assumed to be the mesh file
    int baseId = 1, cellDim, physDim;

    cg_base_read(fileId, baseId, name, &cellDim, &physDim);

    if(cellDim != 2)
        throw Exception("CgnsUnstructuredGrid", "CgnsUnstructuredGrid", "cell dimension must be 2.");

    printf("Reading mesh base \"%s\"...\n", name);

    //- Only a single zone is currently supported. This may change later. If multiple zones, only the first one is read
    int zoneId = 1;

    //- Check zone type, should be unstructured
    ZoneType_t zoneType;
    cg_zone_type(fileId, baseId, zoneId, &zoneType);

    if(zoneType != Unstructured)
        throw Exception("CgnsUnstructuredGrid", "CgnsUnstructuredGrid", "zone type must be unstructured.");

    cgsize_t sizes[2];
    cg_zone_read(fileId, baseId, zoneId, name, sizes);
    nodes_.reserve(sizes[0]);
    cells_.reserve(sizes[1]);

    printf("Loading zone \"%s\" with %d nodes and %d cells...\n", name, sizes[0], sizes[1]);

    //- Read all zone relevant data
    readNodes(fileId, baseId, zoneId, sizes[0], convertToMeters);
    readElements(fileId, baseId, zoneId);
    readBoundaries(fileId, baseId, zoneId);

    cg_close(fileId);

    initConnectivity();
    computeBoundingBox();
}

void CgnsUnstructuredGrid::newFile(const std::string &filename, const std::string &baseName)
{
    if(fileIsOpen_)
        closeFile();

    cg_open(filename.c_str(), CG_MODE_WRITE, &fileId_);
    cg_base_write(fileId_, baseName.c_str(), 2, 2, &baseId_);
    fileIsOpen_ = true;
}

void CgnsUnstructuredGrid::openFile(const std::string &filename)
{
    if(fileIsOpen_)
        closeFile();

    cg_open(filename.c_str(), CG_MODE_MODIFY, &fileId_);

    int nBases;
    cg_nbases(fileId_, &nBases);

    if(nBases != 1)
        throw Exception("CgnsUnstructuredGrid", "openFile", "only one base is allowed.");

    baseId_ = 1;
    fileIsOpen_ = true;
}

void CgnsUnstructuredGrid::closeFile()
{
    cg_close(fileId_);
    fileIsOpen_ = false;
}

int CgnsUnstructuredGrid::addZone(const std::string &zoneName, int nNodes, int nCells)
{
    if(!fileIsOpen_)
        throw Exception("CgnsUnstructuredGrid", "addZone", "no cgns file is currently open.");

    cgsize_t sizes[] = {nNodes, nCells, 0};
    int zoneId;
    cg_zone_write(fileId_, baseId_, zoneName.c_str(), sizes, Unstructured, &zoneId);

    return zoneId;
}

void CgnsUnstructuredGrid::addNodes(int zoneId, const std::vector<Point2D> &nodes)
{
    if(!fileIsOpen_)
        throw Exception("CgnsUnstructuredGrid", "addNodes", "no cgns file is currently open.");

    std::vector<double> coordsX, coordsY;
    coordsX.reserve(nodes.size());
    coordsY.reserve(nodes.size());

    for(const Point2D& node: nodes)
    {
        coordsX.push_back(node.x);
        coordsY.push_back(node.y);
    }

    int coordId;
    cg_coord_write(fileId_, baseId_, zoneId, RealDouble, "CoordinateX", coordsX.data(), &coordId);
    cg_coord_write(fileId_, baseId_, zoneId, RealDouble, "CoordinateY", coordsY.data(), &coordId);
}

int CgnsUnstructuredGrid::addTriCells(int zoneId, const std::vector<cgsize_t> &cells)
{
    if(!fileIsOpen_)
        throw Exception("CgnsUnstructuredGrid", "addTriCells", "no cgns file currently open.");

    int secId;
    cg_section_write(fileId_, baseId_, zoneId, "GridElements", TRI_3, 1, cells.size()/3, 0, cells.data(), &secId);
    return secId;
}

int CgnsUnstructuredGrid::addMixedCells(int zoneId, int nCells, const std::vector<cgsize_t> &cells)
{
    if(!fileIsOpen_)
        throw Exception("CgnsUnstructuredGrid", "addMixedCells", "no cgns file is currently open.");

    int secId;

    cg_section_write(fileId_, baseId_, zoneId, "GridElements", MIXED, 1, nCells, 0, cells.data(), &secId);
    return secId;
}

int CgnsUnstructuredGrid::addBc(int zoneId, const std::string &name, const std::vector<cgsize_t> &faces)
{
    if(!fileIsOpen_)
        throw Exception("CgnsUnstructuredGrid", "addBc", "no cgns file is currently open.");

    int secId;

    int nSecs;
    cg_nsections(fileId_, baseId_, zoneId, &nSecs);

    cgsize_t maxElement = 0, start, end;
    for(int secNo = 1; secNo <= nSecs; ++secNo)
    {
        char secName[256];
        ElementType_t type;
        int nBoundary;

        cg_section_read(fileId_, baseId_, zoneId, secNo, secName, &type, &start, &end, &nBoundary, NULL);
        maxElement = std::max(maxElement, end);
    }

    start = maxElement + 1;
    end = start + faces.size()/2 - 1;

    cg_section_write(fileId_, baseId_, zoneId, (name + "Faces").c_str(), BAR_2, start, end, 0, faces.data(), &secId);

    int bcId;
    cgsize_t pointRange[] = {start, end};

    cg_boco_write(fileId_, baseId_, zoneId, name.c_str(), BCGeneral, PointRange, end - start + 1, pointRange, &bcId);

    return bcId;
}

int CgnsUnstructuredGrid::connectZones(int zoneId, const std::vector<cgsize_t> &faces, int donorZoneId, const std::vector<cgsize_t>& donorCells)
{
    using namespace std;

    if(!fileIsOpen_)
        throw Exception("CgnsUnstructuredGrid", "connectZones", "no cgns file is currently open.");

    int interfaceId;
    cgsize_t sizes[3];
    char donorName[256];

    cg_zone_read(fileId_, baseId_, donorZoneId, donorName, sizes);

    cg_conn_write(fileId_, baseId_, zoneId,
                  (string(donorName) + "_interface").c_str(), EdgeCenter, Abutting1to1,
                  PointList, faces.size(), faces.data(),
                  donorName, Unstructured, CellListDonor,
                  LongInteger, donorCells.size(), donorCells.data(), &interfaceId);

    return interfaceId;
}

//- Private helper methods

void CgnsUnstructuredGrid::readNodes(int fileId, int baseId, int zoneId, int nNodes, Scalar convertToMeters)
{
    std::vector<double> xCoords(nNodes), yCoords(nNodes);
    cgsize_t rmin = 1, rmax = nNodes;

    cg_coord_read(fileId, baseId, zoneId, "CoordinateX", RealDouble, &rmin, &rmax, xCoords.data());
    cg_coord_read(fileId, baseId, zoneId, "CoordinateY", RealDouble, &rmin, &rmax, yCoords.data());

    for(int i = 0; i < nNodes; ++i)
        addNode(Point2D(xCoords[i]*convertToMeters, yCoords[i]*convertToMeters));

    //initNodes();
}

void CgnsUnstructuredGrid::readElements(int fileId, int baseId, int zoneId)
{
    int nSections;
    cg_nsections(fileId, baseId, zoneId, &nSections);

    for(int secId = 1; secId <= nSections; ++secId)
    {
        char name[256];
        ElementType_t type;
        int start, end, nBoundary, parentFlag;

        cg_section_read(fileId, baseId, zoneId, secId, name, &type, &start, &end, &nBoundary, &parentFlag);

        const cgsize_t nElems = end - start + 1;
        cgsize_t nElemNodes;

        switch(type)
        {
        case TRI_3:
            printf("Initializing triangular elements of section \"%s\"...\n", name);
            nElemNodes = 3;
            break;

        case QUAD_4:
            printf("Initializing quadrilateral elements of section \"%s\"...\n", name);
            nElemNodes = 4;
            break;

        case BAR_2:
            continue;

        default:
            throw Exception("CgnsUnstructuredGrid", "readElements", "unsupported element type. Only TRI_3, QUAD_4 and BAR_2 are currently valid.");
        };

        std::vector<cgsize_t> elems(nElemNodes*nElems);
        cg_elements_read(fileId, baseId, zoneId, secId, elems.data(), NULL);

        for(int i = 0; i < nElems; ++i)
        {
            std::vector<Label> nodeIds;

            for(int j = 0; j < nElemNodes; ++j)
                nodeIds.push_back(elems[nElemNodes*i + j] - 1);

            createCell(nodeIds);
        }
    } // end for

    //initCells();
}

void CgnsUnstructuredGrid::readBoundaries(int fileId, int baseId, int zoneId)
{
    using namespace std;

    int nSections;
    cg_nsections(fileId, baseId, zoneId, &nSections);

    map<pair<int, int>, vector<cgsize_t> > eleMap;

    //- Find the boundary elements
    for(int secId = 1; secId <= nSections; ++secId)
    {
        char name[256];
        ElementType_t type;
        int start, end, nBoundary, parentFlag;

        cg_section_read(fileId, baseId, zoneId, secId, name, &type, &start, &end, &nBoundary, &parentFlag);

        switch(type)
        {
        case TRI_3: case QUAD_4:
            continue;

        case BAR_2:
            break;
        }

        const cgsize_t nElems = end - start + 1;
        vector<cgsize_t> elems(2*nElems);

        cg_elements_read(fileId, baseId, zoneId, secId, elems.data(), NULL);
        eleMap.insert(make_pair(make_pair(start, end), elems));
    }

    //- Read the boundaries
    int nBcs;
    cg_nbocos(fileId, baseId, zoneId, &nBcs);

    for(int bcId = 1; bcId <= nBcs; ++bcId)
    {
        char name[256];
        BCType_t bcType;
        PointSetType_t pointSetType;
        cgsize_t nElems;
        int normalIndex;
        cgsize_t normalListSize;
        DataType_t dataType;
        int nDataSet;

        cg_boco_info(fileId, baseId, zoneId, bcId, name, &bcType, &pointSetType, &nElems, &normalIndex, &normalListSize, &dataType, &nDataSet);

        printf("\nCreating boundary patch \"%s\"...\n", name);
        printf("BC type: %s\n", BCTypeName[bcType]);
        printf("Number of faces: %d\n", (int)nElems);

        vector<cgsize_t> elemIds(nElems);
        cg_boco_read(fileId, baseId, zoneId, bcId, elemIds.data(), NULL);

        vector<Label> faces;
        for(cgsize_t elemId: elemIds)
        {
            Label n1, n2;

            for(const auto &entry: eleMap)
            {
                const auto &range = entry.first;

                if(elemId >= range.first && elemId <= range.second)
                {
                    const auto &elems = entry.second;
                    const int i = elemId - range.first;

                    n1 = elems[2*i] - 1;
                    n2 = elems[2*i + 1] - 1;
                    break;
                }
            }

            faces.push_back(findFace(n1, n2));
        }

        applyPatch(name, faces);
    }
}

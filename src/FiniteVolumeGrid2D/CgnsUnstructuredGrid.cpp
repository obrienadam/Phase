#include <cgnslib.h>

#include "CgnsUnstructuredGrid.h"
#include "Exception.h"

CgnsUnstructuredGrid::CgnsUnstructuredGrid(const Input &input)
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

        vector< Ref<Face> > faces;
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

            faces.push_back(ref(faces_[findFace(n1, n2)]));
        }

        applyPatch(name, faces);
    }
}

#include "CgnsUnstructuredGrid.h"
#include "Communicator.h"
#include "Exception.h"

CgnsUnstructuredGrid::CgnsUnstructuredGrid()
        :
        FiniteVolumeGrid2D()
{

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

    if (cellDim != 2)
        throw Exception("CgnsUnstructuredGrid", "CgnsUnstructuredGrid", "cell dimension must be 2.");

    printf("Reading mesh base \"%s\"...\n", name);

    //- Only a single zone is currently supported. This may change later. If multiple zones, only the first one is read
    int zoneId = 1;

    //- Check zone type, should be unstructured
    CGNS_ENUMT(ZoneType_t) zoneType;
    cg_zone_type(fileId, baseId, zoneId, &zoneType);

    if (zoneType != CGNS_ENUMV(Unstructured))
        throw Exception("CgnsUnstructuredGrid", "CgnsUnstructuredGrid", "zone type must be unstructured.");

    cgsize_t sizes[2];
    cg_zone_read(fileId, baseId, zoneId, name, sizes);
    nodes_.reserve(sizes[0]);
    cells_.reserve(sizes[1]);

    printf("Loading zone \"%s\" with %d nodes and %d cells...\n", name, sizes[0], sizes[1]);

    //- Read all zone relevant data
    readNodes(fileId, baseId, zoneId, sizes[0], convertToMeters, input.caseInput().get<std::string>("Grid.origin", "(0,0)"));
    readElements(fileId, baseId, zoneId);
    readBoundaries(fileId, baseId, zoneId);

    cg_close(fileId);

    initConnectivity();
    computeBoundingBox();
}

void CgnsUnstructuredGrid::loadPartitionedGrid(std::shared_ptr<Communicator> comm)
{
    comm_ = comm;

    int fid;
    cg_open(("solution/Proc" + std::to_string(comm_->rank()) + "/Grid.cgns").c_str(), CG_MODE_READ, &fid);

    int bid = 1, zid = 1, dim[2];
    char name[256];
    cg_base_read(fid, bid, name, &dim[0], &dim[1]);

    CGNS_ENUMT(ZoneType_t) zoneType;
    cg_zone_type(fid, bid, zid, &zoneType);
    cgsize_t sizes[2];
    cg_zone_read(fid, bid, zid, name, sizes);

    nodes_.reserve(sizes[0]);
    cells_.reserve(sizes[1]);

    readNodes(fid, bid, zid, sizes[0], 1., Point2D(0., 0.));
    readElements(fid, bid, zid);
    readBoundaries(fid, bid, zid);
    initConnectivity();
    computeBoundingBox();

    std::vector<int> procNo(cells_.size()), globalIds(cells_.size());
    cgsize_t rmin = 1, rmax = cells_.size();

    cg_field_read(fid, bid, zid, 1, "ProcNo", CGNS_ENUMV(Integer), &rmin, &rmax, procNo.data());
    cg_field_read(fid, bid, zid, 1, "GlobalID", CGNS_ENUMV(Integer), &rmin, &rmax, globalIds.data());

    cg_close(fid);

    //- Construct the buffer zones
    sendCellGroups_.resize(comm_->nProcs());
    bufferCellZones_.resize(comm_->nProcs());

    for (int proc = 0; proc < comm_->nProcs(); ++proc)
    {
        sendCellGroups_[proc] = CellGroup("Proc" + std::to_string(proc));
        bufferCellZones_[proc] = CellZone("Proc" + std::to_string(proc), localActiveCells_.registry());
    }

    //- Construct the buffer regions
    for (const Cell &cell: cells_)
        if (procNo[cell.id()] != comm_->rank())
            bufferCellZones_[procNo[cell.id()]].add(cell);

    //- Create a global to local id map
    std::unordered_map<int, int> globalToLocalIdMap;
    for (int id = 0; id < globalIds.size(); ++id)
        globalToLocalIdMap[globalIds[id]] = id;

    //- Communicate send orders
    std::vector<std::vector<unsigned long>> recvOrders(comm_->nProcs());
    for (int proc = 0; proc < comm_->nProcs(); ++proc)
    {
        std::transform(bufferCellZones_[proc].begin(), bufferCellZones_[proc].end(),
                       std::back_inserter(recvOrders[proc]), [&globalIds](const Cell &cell) {
                    return globalIds[cell.id()];
                });

        comm_->isend(proc, recvOrders[proc], proc);
    }

    for (int proc = 0; proc < comm_->nProcs(); ++proc)
    {
        std::vector<unsigned long> sendOrder(comm_->probeSize<unsigned long>(proc, comm_->rank()));
        comm_->recv(proc, sendOrder, comm_->rank());

        for (Label gid: sendOrder)
            sendCellGroups_[proc].add(cells_[globalToLocalIdMap[gid]]);
    }

    comm_->waitAll();
    computeGlobalOrdering();
}

//- Private helper methods

void CgnsUnstructuredGrid::readNodes(int fileId,
                                     int baseId,
                                     int zoneId,
                                     int nNodes,
                                     Scalar convertToMeters,
                                     const Point2D &origin)
{
    std::vector<double> xCoords(nNodes), yCoords(nNodes);
    cgsize_t rmin = 1, rmax = nNodes;

    cg_coord_read(fileId, baseId, zoneId, "CoordinateX", CGNS_ENUMV(RealDouble), &rmin, &rmax, xCoords.data());
    cg_coord_read(fileId, baseId, zoneId, "CoordinateY", CGNS_ENUMV(RealDouble), &rmin, &rmax, yCoords.data());

    for (int i = 0; i < nNodes; ++i)
        addNode((Point2D(xCoords[i], yCoords[i]) + origin) * convertToMeters);
}

void CgnsUnstructuredGrid::readElements(int fileId, int baseId, int zoneId)
{
    int nSections;
    cg_nsections(fileId, baseId, zoneId, &nSections);

    for (int secId = 1; secId <= nSections; ++secId)
    {
        char name[256];
        CGNS_ENUMT(ElementType_t) type;
        int start, end, nBoundary, parentFlag;

        cg_section_read(fileId, baseId, zoneId, secId, name, &type, &start, &end, &nBoundary, &parentFlag);

        const cgsize_t nElems = end - start + 1;
        cgsize_t nElemNodes;

        switch (type)
        {
            case CGNS_ENUMV(TRI_3):
                comm_->printf("Initializing triangular elements of section \"%s\"...\n", name);
                nElemNodes = 3;
                break;

            case CGNS_ENUMV(QUAD_4):
                comm_->printf("Initializing quadrilateral elements of section \"%s\"...\n", name);
                nElemNodes = 4;
                break;

            case CGNS_ENUMV(MIXED): //- Special behaviour needed
                comm_->printf("Initializing mixed elements of section \"%s\"...\n", name);
                {
                    cgsize_t size;
                    cg_ElementDataSize(fileId, baseId, zoneId, secId, &size);
                    std::vector<cgsize_t> elems(size);
                    cg_elements_read(fileId, baseId, zoneId, secId, elems.data(), NULL);

                    for (int i = 0; i < size;)
                    {
                        int n;
                        switch (elems[i++])
                        {
                            case CGNS_ENUMV(TRI_3):
                                n = 3;
                                break;
                            case CGNS_ENUMV(QUAD_4):
                                n = 4;
                                break;
                            default:
                                throw Exception("CgnsUnstructuredGrid", "readElements",
                                                "unsupported mixed element type. Only TRI_3 and QUAD_4 are currently valid.");
                        }

                        std::vector<Label> nodeIds;

                        for (int j = 0; j < n; ++j)
                            nodeIds.push_back(elems[i++] - 1);
                        createCell(nodeIds);
                    }
                }
                continue;

            case CGNS_ENUMV(BAR_2):
                continue;

            default:
                throw Exception("CgnsUnstructuredGrid", "readElements",
                                "unsupported element type. Only TRI_3, QUAD_4, MIXED and BAR_2 are currently valid.");
        };

        std::vector<cgsize_t> elems(nElemNodes * nElems);
        cg_elements_read(fileId, baseId, zoneId, secId, elems.data(), NULL);

        for (int i = 0; i < nElems; ++i)
        {
            std::vector<Label> nodeIds;

            for (int j = 0; j < nElemNodes; ++j)
                nodeIds.push_back(elems[nElemNodes * i + j] - 1);

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
    for (int secId = 1; secId <= nSections; ++secId)
    {
        char name[256];
        CGNS_ENUMT(ElementType_t) type;
        int start, end, nBoundary, parentFlag;

        cg_section_read(fileId, baseId, zoneId, secId, name, &type, &start, &end, &nBoundary, &parentFlag);

        switch (type)
        {
            case CGNS_ENUMV(TRI_3):
            case CGNS_ENUMV(QUAD_4):
            case CGNS_ENUMV(MIXED):
                continue;

            case CGNS_ENUMV(BAR_2):
                break;
        }

        const cgsize_t nElems = end - start + 1;
        vector<cgsize_t> elems(2 * nElems);

        cg_elements_read(fileId, baseId, zoneId, secId, elems.data(), NULL);
        eleMap.insert(make_pair(make_pair(start, end), elems));
    }

    //- Read the boundaries
    int nBcs;
    cg_nbocos(fileId, baseId, zoneId, &nBcs);

    for (int bcId = 1; bcId <= nBcs; ++bcId)
    {
        char name[256];
        CGNS_ENUMT(BCType_t) bcType;
        CGNS_ENUMT(PointSetType_t) pointSetType;
        cgsize_t nElems;
        int normalIndex;
        cgsize_t normalListSize;
        CGNS_ENUMT(DataType_t) dataType;
        int nDataSet;

        cg_boco_info(fileId, baseId, zoneId, bcId, name, &bcType, &pointSetType, &nElems, &normalIndex, &normalListSize,
                     &dataType, &nDataSet);

        printf("\nCreating boundary patch \"%s\"...\n", name);
        printf("BC type: %s\n", BCTypeName[bcType]);
        printf("Number of faces: %d\n", (int) nElems);

        vector<cgsize_t> elemIds(nElems);
        cg_boco_read(fileId, baseId, zoneId, bcId, elemIds.data(), NULL);

        vector<Label> faces;
        for (cgsize_t elemId: elemIds)
        {
            Label n1, n2;

            for (const auto &entry: eleMap)
            {
                const auto &range = entry.first;

                if (elemId >= range.first && elemId <= range.second)
                {
                    const auto &elems = entry.second;
                    const int i = elemId - range.first;

                    n1 = elems[2 * i] - 1;
                    n2 = elems[2 * i + 1] - 1;
                    break;
                }
            }

            faces.push_back(findFace(n1, n2));
        }

        createPatch(name, faces);
    }
}

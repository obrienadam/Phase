#include <cgnslib.h>

#include "CgnsUnstructuredQuadGrid.h"
#include "Exception.h"

CgnsUnstructuredQuadGrid::CgnsUnstructuredQuadGrid(const Input &input)
{
    std::string filename = input.caseInput().get<std::string>("Grid.filename");

    int fileId;
    char name[256];

    cg_open(filename.c_str(), CG_MODE_READ, &fileId);

    //- The first base is assumed to be the mesh file
    int baseId = 1, cellDim, physDim;

    cg_base_read(fileId, baseId, name, &cellDim, &physDim);

    if(cellDim != 2)
        throw Exception("CgnsUnstructuredQuadGrid", "CgnsUnstructuredQuadGrid", "cell dimension must be 2.");

    printf("Reading mesh base \"%s\"...\n", name);

    //- Only a single zone is currently supported. This may change later. If multiple zones, only the first one is read
    int zoneId = 1;

    //- Check zone type, should be unstructured
    ZoneType_t zoneType;
    cg_zone_type(fileId, baseId, zoneId, &zoneType);

    if(zoneType != Unstructured)
        throw Exception("CgnsUnstructuredQuadGrid", "CgnsUnstructuredQuadGrid", "zone type must be unstructured.");

    cgsize_t sizes[2];

    cg_zone_read(fileId, baseId, zoneId, name, sizes);

    printf("Loading zone \"%s\" with %d nodes and %d cells...\n", name, sizes[0], sizes[1]);

    nodes_.reserve(sizes[0]);
    cells_.reserve(sizes[1]);

    //- Read the node coordinates
    std::vector<Scalar> coordsX(sizes[0]), coordsY(sizes[0]);
    int one = 1;

    cg_coord_read(fileId, baseId, zoneId, "CoordinateX", RealDouble, &one, &sizes[0], coordsX.data());
    cg_coord_read(fileId, baseId, zoneId, "CoordinateY", RealDouble, &one, &sizes[0], coordsY.data());

    //- Read the cell connectivity
    int nSections;
    cg_nsections(fileId, baseId, zoneId, &nSections);

    std::vector<cgsize_t> quadElements;
    std::pair<cgsize_t, cgsize_t> quadElementRange;

    std::vector<cgsize_t> triElements;
    std::pair<cgsize_t, cgsize_t> triElementRange;

    std::vector<cgsize_t> barElements;
    std::pair<cgsize_t, cgsize_t> barElementRange;

    for(int sec = 1; sec <= nSections; ++sec)
    {
        ElementType_t type;
        int nBoundary, parentFlag;
        cgsize_t eBeg, eEnd;

        cg_section_read(fileId, baseId, zoneId, sec, name, &type, &eBeg, &eEnd, &nBoundary, &parentFlag);

        switch(type)
        {
        case QUAD_4:
            quadElementRange.first = eBeg - 1; // converted to C array notation , range = [first, second)
            quadElementRange.second = eEnd;
            quadElements.resize(4*(quadElementRange.second - quadElementRange.first));

            cg_elements_read(fileId, baseId, zoneId, sec, quadElements.data(), NULL);
            break;

        case TRI_3:
            triElementRange.first = eBeg - 1;
            triElementRange.second = eEnd;
            triElements.resize(3*(triElementRange.second - triElementRange.first));

            cg_elements_read(fileId, baseId, zoneId, sec, triElements.data(), NULL);
            break;

            break;

        case BAR_2:
            barElementRange.first = eBeg - 1;
            barElementRange.second = eEnd;
            barElements.resize(2*(barElementRange.second - barElementRange.first));

            cg_elements_read(fileId, baseId, zoneId, sec, barElements.data(), NULL);
            break;

        default:
            throw Exception("CgnsUnstructuredGrid", "CgnsUnstructuredGrid", "Element types must either be \"QUAD_4\", \"TRI_3\" or \"BAR_2\".");
        };
    }

    //- Initialize nodes and cells (must be done before patches are constructed)
    for(int i = 0; i < sizes[0]; ++i)
        addNode(Point2D(coordsX[i], coordsY[i]));

    for(int i = 0, nQuads = quadElements.size()/4; i < nQuads; ++i)
    {
        std::vector<Label> nodeIds;
        nodeIds.push_back(quadElements[4*i] - 1);
        nodeIds.push_back(quadElements[4*i + 1] - 1);
        nodeIds.push_back(quadElements[4*i + 2] - 1);
        nodeIds.push_back(quadElements[4*i + 3] - 1);

        createCell(nodeIds);
    }

    for(int i = 0, nTris = triElements.size()/3; i < nTris; ++i)
    {
        std::vector<Label> nodeIds;
        nodeIds.push_back(triElements[3*i] - 1);
        nodeIds.push_back(triElements[3*i + 1] - 1);
        nodeIds.push_back(triElements[3*i + 2] - 1);

        createCell(nodeIds);
    }

    //- Initialize boundary patches
    int nBcs;
    cg_nbocos(fileId, baseId, zoneId, &nBcs);

    for(int bd = 1; bd <= nBcs; ++bd)
    {
        GridLocation_t bcLoc;
        cg_boco_gridlocation_read(fileId, baseId, zoneId, bd, &bcLoc);

        if(bcLoc != EdgeCenter)
            throw Exception("CgnsUnstructuredQuadGrid", "CgnsUnstructuredQuadGrid", "all boundary conditions must be edge-centered.");

        //- Get the bc info
        BCType_t bcType;
        PointSetType_t ptSetType;
        cgsize_t nPts, normalListSize;
        int normalIndex, nDataSet;
        DataType_t normalDataType;

        cg_boco_info(fileId,  baseId, zoneId, bd, name, &bcType, &ptSetType, &nPts, &normalIndex, &normalListSize, &normalDataType, &nDataSet);

        printf("\nCreating boundary patch \"%s\"...\n", name);
        printf("BC type: %s\n", BCTypeName[bcType]);
        printf("N faces: %d\n", (int)nPts);

        std::vector<cgsize_t> boundaryElements(nPts);

        cg_boco_read(fileId, baseId, zoneId, bd, boundaryElements.data(), NULL);

        std::vector< Ref<Face> > faces;
        for(cgsize_t ele: boundaryElements)
        {
            Label id = ele - 1 - barElementRange.first;
            Label n1 = barElements[2*id] - 1, n2 = barElements[2*id + 1] - 1;

            Face& face = faces_[findFace(n1, n2)];

            faces.push_back(std::ref(face));
        }

        applyPatch(name, faces);
    }

    cg_close(fileId);
    initNodes();
    initCells();

    //- Too a quick validity check!
    for(const Face& face: faces_)
    {
        if(face.isBoundary() && !face.belongsToPatch())
            throw Exception("CgnsUnstructuredGrid", "CgnsUnstructuredGrid", "One or more boundary faces has not been associated with a boundary patch.");
    }

    computeBoundingBox();
}

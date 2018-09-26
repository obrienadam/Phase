#include "System/CgnsFile.h"

#include "CgnsUnstructuredGrid.h"

CgnsUnstructuredGrid::CgnsUnstructuredGrid()
        :
        FiniteVolumeGrid2D()
{

}

CgnsUnstructuredGrid::CgnsUnstructuredGrid(const Input &input)
        :
        CgnsUnstructuredGrid()
{
    load(input.caseInput().get<std::string>("Grid.filename"),
         input.caseInput().get<std::string>("Grid.origin", "(0,0)"));
}

void CgnsUnstructuredGrid::load(const std::string &filename, const Point2D &origin)
{
    CgnsFile file(filename, CgnsFile::READ);

    auto base = file.readBase(1);

    //- Check to make sure it is a workable base
    comm_->printf("Read base \"%s\".\n", base.name.c_str());

    if (base.cellDim != 2)
        throw Exception("CgnsUnstructuredGrid", "CgnsUnstructuredGrid", "cell dimension must be be 2.");

    auto zone = file.readZone(1, 1);

    if (zone.type != "Unstructured")
        throw Exception("CgnsUnstructuredGrid", "CgnsUnstructuredGrid", "zone type must be \"Unstructured\".");

    int nSections = file.nSections(1, 1);
    std::vector<std::vector<int>> adjncy;

    //- Gather all elements, including the boundary patch elements
    for (int sid = 1; sid <= nSections; ++sid)
    {
        auto section = file.readSection(1, 1, sid);

        if (adjncy.size() < section.end)
            adjncy.resize(section.end);
        for (auto id = section.start; id <= section.end; ++id)
            adjncy[id - 1].assign(section.cind.begin() + section.cptr[id - section.start],
                                  section.cind.begin() + section.cptr[id - section.start + 1]);
    }

    std::vector<int> eptr(1, 0), eind;

    for (const auto &elem: adjncy)
    {
        eptr.push_back(eptr.back() + elem.size());
        eind.insert(eind.end(), elem.begin(), elem.end());
    }

    std::transform(eind.begin(), eind.end(), eind.begin(), [](int id)
    { return id - 1; });

    std::vector<Point2D> nodes = file.readCoords<Point2D>(1, 1);
    std::transform(nodes.begin(), nodes.end(), nodes.begin(), [origin](const Point2D &node)
    { return node + origin; });

    std::unordered_map<std::string, std::vector<Label>> patches;
    int nBoCos = file.nBoCos(1, 1);

    for (int bcid = 1; bcid <= nBoCos; ++bcid)
    {
        auto boco = file.readBoCo(1, 1, bcid);
        std::vector<Label> pts;

        if (boco.pointListType == "PointRange")
        {
            pts.insert(pts.end(),
                       eind.begin() + eptr[boco.pnts[0] - 1],
                       eind.begin() + eptr[boco.pnts[1] - 1 + 1]);
        }
        else if (boco.pointListType == "PointList")
        {
            for (auto pnt: boco.pnts)
            {
                pts.insert(pts.end(),
                           eind.begin() + eptr[pnt - 1],
                           eind.begin() + eptr[pnt - 1 + 1]);
            }
        }
        else
        {
            throw Exception("CgnsUnstructuredGrid", "CgnsUnstructuredGrid",
                            "bad point set type \"" + boco.pointListType + "\"");
        }

        patches[boco.name] = std::move(pts);

        comm_->printf("Read boundary patch \"%s\".\n", boco.name.c_str());
    }

    file.close();

    //- init all other elements
    std::vector<Label> cptr(1, 0), cind;

    for (auto id = 0; id < eptr.size() - 1; ++id)
        if (eptr[id + 1] - eptr[id] > 2)
        {
            cptr.push_back(cptr.back() + eptr[id + 1] - eptr[id]);
            cind.insert(cind.end(), eind.begin() + eptr[id], eind.begin() + eptr[id + 1]);
        }

    init(nodes, cptr, cind, Point2D(0., 0.));
    initPatches(patches);
}

void CgnsUnstructuredGrid::readPartitionData(const std::string &filename)
{
    CgnsFile file(filename);

    int nSols = file.nSolutions(1, 1);

    for (int sid = 1; sid <= nSols; ++sid)
    {
        auto soln = file.readSolution(1, 1, 1);

        if (soln.name == "Info" || soln.name == "info")
        {
            auto procNo = file.readField<int>(1, 1, sid, 1, cells_.size(), "ProcNo").data;
            auto globalIds = file.readField<int>(1, 1, sid, 1, cells_.size(), "GlobalID").data;

            initCommBuffers(std::vector<Label>(procNo.begin(), procNo.end()),
                            std::vector<Label>(globalIds.begin(), globalIds.end()));

            break;
        }
    }

    file.close();
}

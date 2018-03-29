#include <numeric>

#include <cgnslib.h>

#include "CgnsFile.h"

CgnsFile::CgnsFile(const std::string &filename, Mode mode)
{
    open(filename, mode);
}

CgnsFile::~CgnsFile()
{
    close();
}

int CgnsFile::open(const std::string &filename, Mode mode)
{
    switch (mode)
    {
        case READ:
            cg_open(filename.c_str(), CG_MODE_READ, &_fid);
            break;

        case WRITE:
            cg_open(filename.c_str(), CG_MODE_WRITE, &_fid);
            break;

        case MODIFY:
            cg_open(filename.c_str(), CG_MODE_MODIFY, &_fid);
            break;
    }

    return _fid;
}

void CgnsFile::close()
{
    cg_close(_fid);
}

int CgnsFile::createBase(const std::string &basename, int cellDim, int physDim)
{
    int bid;
    cg_base_write(_fid, basename.c_str(), cellDim, physDim, &bid);
    return bid;
}

std::tuple<int, int> CgnsFile::writeCoordinates(int bid, int zid, const std::vector<Point2D> &coords)
{
    std::vector<double> x(coords.size()), y(coords.size());
    std::transform(coords.begin(), coords.end(), x.begin(), [](const Point2D &pt)
    { return pt.x; });
    std::transform(coords.begin(), coords.end(), y.begin(), [](const Point2D &pt)
    { return pt.y; });

    std::tuple<int, int> cid;
    cg_coord_write(_fid, bid, zid, CGNS_ENUMV(RealDouble), "CoordinateX", x.data(), &std::get<0>(cid));
    cg_coord_write(_fid, bid, zid, CGNS_ENUMV(RealDouble), "CoordinateY", y.data(), &std::get<1>(cid));
    return cid;
}

int CgnsFile::createStructuredZone(int bid, const std::string &zonename,
                                   int nNodesI, int nNodesJ,
                                   int nCellsI, int nCellsJ)
{
    int zid;
    cgsize_t size[] = {
            nNodesI, nNodesJ, nCellsI, nCellsJ, 0, 0
    };

    cg_zone_write(_fid, bid, zonename.c_str(), size, CGNS_ENUMV(Structured), &zid);

    return zid;
}

int CgnsFile::createStructuredZone(int bid,
                                   const std::string &zonename,
                                   int nNodesI, int nNodesJ, int nNodesK,
                                   int nCellsI, int nCellsJ, int nCellsK)
{
    int zid;
    cgsize_t size[] = {
            nNodesI, nNodesJ, nNodesK, nCellsI, nCellsJ, nCellsK, 0, 0, 0
    };

    cg_zone_write(_fid, bid, zonename.c_str(), size, CGNS_ENUMV(Structured), &zid);

    return zid;
}

int CgnsFile::createUnstructuredZone(int bid, const std::string &zonename, int nNodes, int nCells)
{
    int zid;
    cgsize_t size[] = {nNodes, nCells, 0};
    cg_zone_write(_fid, bid, zonename.c_str(), size, CGNS_ENUMV(Unstructured), &zid);
    return zid;
}

int CgnsFile::writeMixedElementSection(int bid, int zid, const std::string &sectionname,
                                       int start, int end, const std::pair<std::vector<int>, std::vector<int>> &conn)
{
    std::vector<cgsize_t> elems;
    elems.reserve(5 * (end - start + 1));

    for (auto i = 0; i < conn.first.size() - 1; ++i)
    {
        switch (conn.first[i + 1] - conn.first[i])
        {
            case 2:
                elems.push_back(CGNS_ENUMV(BAR_2));
                break;
            case 3:
                elems.push_back(CGNS_ENUMV(TRI_3));
                break;
            case 4:
                elems.push_back(CGNS_ENUMV(QUAD_4));
                break;
        }

        for (auto j = conn.first[i]; j < conn.first[i + 1]; ++j)
            elems.push_back(conn.second[j] + 1);
    }

    int sid;
    cg_section_write(_fid, bid, zid, sectionname.c_str(), CGNS_ENUMV(MIXED),
                     start, end, 0, elems.data(), &sid);
    return sid;
}

int CgnsFile::writeBarElementSection(int bid, int zid, const std::string &sectionname,
                                     int start, int end, const std::vector<Label> &elements)
{
    std::vector<cgsize_t> elems(elements.size());
    std::transform(elements.begin(), elements.end(), elems.begin(), [](int v)
    { return v + 1; });

    int sid;
    cg_section_write(_fid, bid, zid, sectionname.c_str(), CGNS_ENUMV(BAR_2),
                     start, end, 0, elems.data(), &sid);

    return sid;
}

int CgnsFile::writeQuadElementSection(int bid, int zid, const std::string &sectionname,
                                      int start, int end, std::vector<int> &elements)
{
    std::vector<cgsize_t> elems(elements.size());
    std::transform(elements.begin(), elements.end(), elems.begin(), [](int v)
    { return v + 1; });

    int sid;
    cg_section_write(_fid, bid, zid, sectionname.c_str(), CGNS_ENUMV(QUAD_4),
                     start, end, 0, elems.data(), &sid);
    return sid;
}

int CgnsFile::writeBC(int bid, int zid, const std::string &bcname, int start, int end)
{
    int bcid;
    cgsize_t pnts[] = {start, end};
    cg_boco_write(_fid, bid, zid, bcname.c_str(), CGNS_ENUMV(BCGeneral), CGNS_ENUMV(PointRange), 2, pnts, &bcid);
    cg_boco_gridlocation_write(_fid, bid, zid, bcid, CGNS_ENUMV(EdgeCenter));
    return bcid;
}

void CgnsFile::linkNode(int bid, int zid,
                        const std::string &nodename, const std::string &filename, const std::string &nameInFile)
{
    cg_goto(_fid, bid, "Zone_t", zid, "end");
    cg_link_write(nodename.c_str(), filename.c_str(), nameInFile.c_str());
}

void CgnsFile::linkNode(int bid, int zid, int sid,
                        const std::string &nodename, const std::string &filename, const std::string &nameInFile)
{
    cg_goto(_fid, bid, "Zone_t", zid, "FlowSolution_t", sid, "end");
    cg_link_write(nodename.c_str(), filename.c_str(), nameInFile.c_str());
}

int CgnsFile::writeSolution(int bid, int zid, const std::string &solnname)
{
    int solid;
    cg_sol_write(_fid, bid, zid, solnname.c_str(), CGNS_ENUMV(CellCenter), &solid);
    return solid;
}

int CgnsFile::writeField(int bid, int zid, int sid, const std::string &fieldname, const std::vector<int> &field)
{
    int fid;
    cg_field_write(_fid, bid, zid, sid, CGNS_ENUMV(Integer), fieldname.c_str(), field.data(), &fid);
    return fid;
}

int CgnsFile::writeField(int bid, int zid, int sid, const std::string &fieldname, const std::vector<double> &field)
{
    int fid;
    auto data = field;
    data.push_back(98247389);
    data.push_back(89742384);
    cg_field_write(_fid, bid, zid, sid, CGNS_ENUMV(RealDouble), fieldname.c_str(), data.data(), &fid);
    return fid;
}

int CgnsFile::writeField(int bid, int zid, int sid, const std::string &fieldname, const std::vector<Vector2D> &field)
{
    int fid;
    std::vector<double> xComp(field.size()), yComp(field.size());
    std::transform(field.begin(), field.end(), xComp.begin(), [](const Vector2D &u)
    { return u.x; });
    std::transform(field.begin(), field.end(), yComp.begin(), [](const Vector2D &u)
    { return u.y; });

    writeField(bid, zid, sid, fieldname + "X", xComp);
    writeField(bid, zid, sid, fieldname + "Y", yComp);

    return fid;
}
#include <numeric>

#include <cgnslib.h>

#include "CgnsFile.h"
#include "Exception.h"

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
    close();

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

int CgnsFile::nBases() const
{
    int nbases;
    cg_nbases(_fid, &nbases);
    return nbases;
}

CgnsFile::Base CgnsFile::readBase(int bid) const
{
    Base base;
    char buff[256];

    cg_base_read(_fid, bid, buff, &base.cellDim, &base.physDim);

    base.id = bid;
    base.name = std::string(buff);

    return base;
}

int CgnsFile::nZones(int bid) const
{
    int nzones;
    cg_nzones(_fid, bid, &nzones);
    return nzones;
}

CgnsFile::Zone CgnsFile::readZone(int bid, int zid) const
{
    Zone zone;
    char buff[256];
    CGNS_ENUMT(ZoneType_t) type;
    cg_zone_type(_fid, bid, zid, &type);
    cg_zone_read(_fid, bid, zid, buff, zone.size);

    zone.id = zid;
    zone.name = std::string(buff);
    zone.type = std::string(cg_ZoneTypeName(type));

    return zone;
};

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

std::tuple<int, int, int> CgnsFile::writeCoordinates(int bid, int zid, const std::vector<Point3D> &coords)
{
    std::vector<double> x(coords.size()), y(coords.size()), z(coords.size());

    std::transform(coords.begin(), coords.end(), x.begin(), [](const Point3D &pt)
    { return pt.x; });

    std::transform(coords.begin(), coords.end(), y.begin(), [](const Point3D &pt)
    { return pt.y; });

    std::transform(coords.begin(), coords.end(), z.begin(), [](const Point3D &pt)
    { return pt.z; });

    std::tuple<int, int, int> cid;
    cg_coord_write(_fid, bid, zid, CGNS_ENUMV(RealDouble), "CoordinateX", x.data(), &std::get<0>(cid));
    cg_coord_write(_fid, bid, zid, CGNS_ENUMV(RealDouble), "CoordinateY", y.data(), &std::get<1>(cid));
    cg_coord_write(_fid, bid, zid, CGNS_ENUMV(RealDouble), "CoordinateZ", z.data(), &std::get<2>(cid));
    return cid;
}

template<>
std::vector<Point2D> CgnsFile::readCoords(int bid, int zid) const
{
    char buff[256];
    cgsize_t size[3];
    cg_zone_read(_fid, bid, zid, buff, size);

    std::vector<double> xCoord(size[0]), yCoord(size[0]);

    cgsize_t rmin = 1, rmax = size[0];
    cg_coord_read(_fid, bid, zid, "CoordinateX", CGNS_ENUMV(RealDouble), &rmin, &rmax, xCoord.data());
    cg_coord_read(_fid, bid, zid, "CoordinateY", CGNS_ENUMV(RealDouble), &rmin, &rmax, yCoord.data());

    std::vector<Point2D> coords(size[0]);
    std::transform(xCoord.begin(), xCoord.end(), yCoord.begin(), coords.begin(), [](Scalar x, Scalar y)
    {
        return Point2D(x, y);
    });

    return coords;
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

int CgnsFile::nSections(int bid, int zid) const
{
    int nsections;
    cg_nsections(_fid, bid, zid, &nsections);
    return nsections;
}

CgnsFile::Section CgnsFile::readSection(int bid, int zid, int sid) const
{
    char buff[256];
    CGNS_ENUMT(ElementType_t) type;

    Section section;

    cg_section_read(_fid, bid, zid, sid, buff, &type, &section.start, &section.end, &section.nbndry,
                    &section.parentFlag);

    cgsize_t dataSize;
    cg_ElementDataSize(_fid, bid, zid, sid, &dataSize);

    std::vector<cgsize_t> elements(dataSize);
    cg_elements_read(_fid, bid, zid, sid, elements.data(), nullptr);

    auto getNVerts = [](CGNS_ENUMT(ElementType_t) type)
    {
        switch (type)
        {
            case CGNS_ENUMV(BAR_2):
                return 2;
            case CGNS_ENUMV(TRI_3):
                return 3;
            case CGNS_ENUMV(QUAD_4):
                return 4;
            default:
                throw Exception("CgnsFile",
                                "readSection",
                                "unsupported element type \"" + std::string(cg_ElementTypeName(type)) + "\".\n");
        }
    };

    std::vector<int> eptr(1, 0), eind;

    int nVerts;
    if (type == CGNS_ENUMV(MIXED))
    {
        for (int i = 0; i < elements.size(); i += nVerts)
        {
            nVerts = getNVerts(CGNS_ENUMT(ElementType_t)(elements[i++]));

            eptr.push_back(eptr.back() + nVerts);
            for (int j = 0; j < nVerts; ++j)
                eind.push_back(elements[i + j]);
        }
    }
    else
    {
        nVerts = getNVerts(type);

        for (int i = 0; i < section.end - section.start + 1; ++i)
        {
            eptr.push_back(eptr.back() + nVerts);
            eind.insert(eind.end(), elements.begin() + nVerts * i, elements.begin() + nVerts * (i + 1));
        }
    }

    section.id = sid;
    section.name = std::string(buff);
    section.type = std::string(cg_ElementTypeName(type));
    section.cptr = std::move(eptr);
    section.cind = std::move(eind);

    return section;
}

int CgnsFile::writeMixedElementSection(int bid, int zid, const std::string &sectionname,
                                       int start, int end, const std::vector<int> &eptr, const std::vector<int> &eind)
{
    std::vector<cgsize_t> elements;

    for (auto i = 0; i < eptr.size() - 1; ++i)
    {

        switch (eptr[i + 1] - eptr[i])
        {
            case 2:
                elements.push_back(CGNS_ENUMV(BAR_2));
                break;
            case 3:
                elements.push_back(CGNS_ENUMV(TRI_3));
                break;
            case 4:
                elements.push_back(CGNS_ENUMV(QUAD_4));
                break;
            default:
                throw Exception("CgnsFile", "writeMixedElementSection", "bad element.");
        }

        elements.insert(elements.end(), eind.begin() + eptr[i], eind.begin() + eptr[i + 1]);
    }

    int sid;
    cg_section_write(_fid, bid, zid, sectionname.c_str(), CGNS_ENUMV(MIXED), start, end, 0, elements.data(), &sid);

    return sid;
}

int CgnsFile::writeBarElementSection(int bid, int zid, const std::string &sectionname, int start, int end,
                                     const std::vector<int> &elements)
{
    int sid;
    cg_section_write(_fid, bid, zid, sectionname.c_str(), CGNS_ENUMV(BAR_2), start, end, 0, elements.data(), &sid);
    return sid;
}

int CgnsFile::nBoCos(int bid, int zid) const
{
    int nbcs;
    cg_nbocos(_fid, bid, zid, &nbcs);
    return nbcs;
}

CgnsFile::BoCo CgnsFile::readBoCo(int bid, int zid, int bcid) const
{
    char buff[256];
    CGNS_ENUMT(BCType_t) type;
    CGNS_ENUMT(PointSetType_t) ptSetType;
    int normalIndex;
    cgsize_t normalListSize;
    CGNS_ENUMT(DataType_t) normalDataType;
    int nDataSet;

    cgsize_t npts;

    cg_boco_info(_fid, bid, zid, bcid, buff, &type, &ptSetType, &npts, &normalIndex, &normalListSize, &normalDataType,
                 &nDataSet);

    std::vector<cgsize_t> pnts(npts);
    cg_boco_read(_fid, bid, zid, bcid, pnts.data(), nullptr);

    BoCo boco;

    boco.id = bcid;
    boco.name = std::string(buff);
    boco.type = std::string(cg_BCTypeName(type));
    boco.pointListType = std::string(cg_PointSetTypeName(ptSetType));
    boco.pnts = pnts;

    return boco;
}

int CgnsFile::writeBoCo(int bid, int zid, const std::string &bcname, int start, int end)
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

int CgnsFile::nSolutions(int bid, int zid) const
{
    int nsols;
    cg_nsols(_fid, bid, zid, &nsols);
    return nsols;
}

CgnsFile::Solution CgnsFile::readSolution(int bid, int zid, int sid) const
{
    char buff[256];
    CGNS_ENUMT(GridLocation_t) location;

    Solution soln;

    cg_sol_info(_fid, bid, zid, sid, buff, &location);
    cg_sol_size(_fid, bid, zid, sid, &soln.dataDim, soln.dimVals);

    soln.name = buff;
    soln.location = cg_GridLocationName(location);
    soln.id = sid;

    return soln;
}

int CgnsFile::writeSolution(int bid, int zid, const std::string &solnname)
{
    int solid;
    cg_sol_write(_fid, bid, zid, solnname.c_str(), CGNS_ENUMV(CellCenter), &solid);
    return solid;
}

template<>
CgnsFile::Field<int> CgnsFile::readField(int bid, int zid, int sid, int rmin, int rmax, const std::string &fieldname)
{
    Field<int> field;

    field.name = fieldname;
    field.rmin = {rmin, 1, 1};
    field.rmax = {rmax, 1, 1};
    field.data.resize(rmax - rmin + 1);

    cg_field_read(_fid, bid, zid, sid,
                  field.name.c_str(),
                  CGNS_ENUMV(Integer),
                  field.rmin.data(),
                  field.rmax.data(),
                  field.data.data());

    return field;
}

template<>
CgnsFile::Field<double> CgnsFile::readField(int bid, int zid, int sid, int rmin, int rmax, const std::string &fieldname)
{
    Field<double> field;

    field.name = fieldname;
    field.rmin = {rmin, 1, 1};
    field.rmax = {rmax, 1, 1};
    field.data.resize(rmax - rmin + 1);

    cg_field_read(_fid, bid, zid, sid,
                  field.name.c_str(),
                  CGNS_ENUMV(RealDouble),
                  field.rmin.data(),
                  field.rmax.data(),
                  field.data.data());

    return field;
}

template<>
int CgnsFile::writeField(int bid, int zid, int sid, const std::string &fieldname, const std::vector<int> &field)
{
    int fid;
    cg_field_write(_fid, bid, zid, sid, CGNS_ENUMV(Integer), fieldname.c_str(), field.data(), &fid);
    return fid;
}

template<>
int CgnsFile::writeField(int bid, int zid, int sid, const std::string &fieldname, const std::vector<Label> &field)
{
    return writeField(bid, zid, sid, fieldname, std::vector<int>(field.begin(), field.end()));
}

template<>
int CgnsFile::writeField(int bid, int zid, int sid, const std::string &fieldname, const std::vector<double> &field)
{
    int fid;
    auto data = field;
    cg_field_write(_fid, bid, zid, sid, CGNS_ENUMV(RealDouble), fieldname.c_str(), data.data(), &fid);
    return fid;
}

template<>
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

int CgnsFile::nDescriptorNodes(int bid) const
{
    cg_goto(_fid, bid, "end");
    int ndescriptors;
    cg_ndescriptors(&ndescriptors);
    return ndescriptors;
}

void CgnsFile::writeDescriptorNode(int bid, const std::string &name, const std::string &text)
{
    cg_goto(_fid, bid, "end");
    cg_descriptor_write(name.c_str(), text.c_str());
}

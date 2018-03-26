#include <cgnslib.h>

#include "CgnsFile.h"

void CgnsFile::~CgnsFile()
{
    close();
}

void CgnsFile::open(const std::string &filename, Mode mode)
{
    if (isOpen())
        close();

    switch (mode)
    {
        case READ:
            cg_open(filename.c_str(), CG_MODE_READ, &fid_);
            break;
        case WRITE:
            cg_open(filename.c_str(), CG_MODE_WRITE, &fid_);
            break;
        case MODIFY:
            cg_open(filename.c_str(), CG_MODE_MODIFY, &fid_);
            break;
    }
}

void CgnsFile::close()
{
    if (isOpen())
    {
        cg_close(fid_);
        fid_ = -1;
    }
}

bool CgnsFile::isOpen() const
{
    return fid_ > 0;
}

int CgnsFile::createBase(const std::string &basename, int cellDim, int physDim)
{
    int bid;
    cg_base_write(fid_, basename.c_str(), cellDim, physDim, &bid);
    return bid;
}

int CgnsFile::nBases() const
{
    int nbases;
    cg_nbases(fid_, &nbases);
    return nbases;
}

int CgnsFile::createZone(int bid, const std::string &zonename, size_t nNodes, size_t nCells, ZoneType type)
{
    CG_ZoneType_t zonetype;

    switch (type)
    {
        case STRUCTURED:
            zonetype = CGNS_ENUMT(Structured);
            break;
        case UNSTRUCTURED:
            zonetype = CGNS_ENUMT(Unstructured);
            break;
    }

    cgsize_t sizes[] = {
            nNodes,
            nCells,
            0
    };

    int zid;
    cg_zone_write(fid_, bid, zonename.c_str(), sizes, zonetype, &zid);

    return zid;
}

int CgnsFile::writeCoords(int bid, int zid, const std::vector<Point2D> &coords)
{
    std::vector<Scalar> xCoords(coords.size()), yCoords(coords.size());

    std::transform(coords.begin(), coords.end(), xCoords.begin(), [](const Point2D &pt) { return pt.x; });
    std::transform(coords.begin(), coords.end(), yCoords.begin(), [](const Point2D &pt) { return pt.y; });

    int cid;
    cg_coord_write(fid_, bid, zid, CGNS_ENUMT(RealDouble), "CoordinateX", xCoords.data(), &cid);
    cg_coord_write(fid_, bid, zid, CGNS_ENUMT(RealDouble), "CoordinateY", yCoords.data(), &cid);
}
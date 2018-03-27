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

int CgnsFile::createStructuredZone(int bid, const std::string &zonename,
                                   int nNodesI, int nNodesJ,
                                   int nCellsI, int nCellsJ)
{
    int zid;
    cgsize_t size[] = {
            nNodesI, nNodesJ, nCellsI, nCellsJ, 0, 0
    };

    cg_zone_write(_fid, bid, zonename.c_str(), size, CGNS_ENUMT(Structured), &zid);

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

    cg_zone_write(_fid, bid, zonename.c_str(), size, CGNS_ENUMT(Structured), &zid);

    return zid;
}

#ifndef CGNS_FILE_H
#define CGNS_FILE_H

#include "Point2D.h"

class CgnsFile
{
public:

    enum Mode
    {
        READ, WRITE, MODIFY
    };

    enum ZoneType
    {
        STRUCTURED, UNSTRUCTURED
    };

    ~CgnsFile();

    void open(const std::string &filename, Mode mode);

    void close();

    bool isOpen() const;

    //- Bases
    int createBase(const std::string &basename, int cellDim, int physDim);

    int nBases() const;

    //- Zones
    int createZone(int bid, const std::string &zonename, size_t nNodes, size_t nCells, ZoneType type);

    //- Coordinates
    int writeCoords(int bid, int zid, const std::vector<Point2D> &coords);

protected:

    int fid_;
};

#endif

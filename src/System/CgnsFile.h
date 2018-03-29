#ifndef PHASE_CGNS_FILE_H
#define PHASE_CGNS_FILE_H

#include <string>
#include <vector>
#include <tuple>

#include "2D/Geometry/Point2D.h"

class CgnsFile
{
public:

    enum Mode
    {
        READ, WRITE, MODIFY
    };

    CgnsFile()
    {}

    CgnsFile(const std::string &filename, Mode mode = READ);

    ~CgnsFile();

    int open(const std::string &filename, Mode mode);

    void close();

    int createBase(const std::string &basename, int cellDim, int physDim);

    //- zones
    std::tuple<int, int> writeCoordinates(int bid, int zid, const std::vector<Point2D> &coords);

    int createStructuredZone(int bid, const std::string &zonename,
                             int nNodesI, int nNodesJ,
                             int nCellsI, int nCellsJ);

    int createStructuredZone(int bid, const std::string &zonename,
                             int nNodesI, int nNodesJ, int nNodesK,
                             int nCellsI, int nCellsJ, int nCellsK);

    int createUnstructuredZone(int bid, const std::string &zonename, int nNodes, int nCells);

    int writeMixedElementSection(int bid, int zid, const std::string &sectionname,
                                 int start, int end, const std::pair<std::vector<int>, std::vector<int>> &conn);

    int writeBarElementSection(int bid, int zid, const std::string &sectionname,
                               int start, int end, const std::vector<Label> &elements);

    int writeQuadElementSection(int bid, int zid, const std::string &sectionname,
                                int start, int end, std::vector<int> &elements);

    int writeBC(int bid, int zid, const std::string& bcname, int start, int end);

    //- linking

    void linkNode(int bid, int zid,
                  const std::string &nodename, const std::string &filename, const std::string &nameInFile);

    void linkNode(int bid, int zid, int sid,
                  const std::string& nodename, const std::string& filename, const std::string& nameInFile);

    //- Solutions
    int writeSolution(int bid, int zid, const std::string &solnname);

    //- Fields
    int writeField(int bid, int zid, int sid, const std::string &fieldname, const std::vector<int> &field);

    int writeField(int bid, int zid, int sid, const std::string &fieldname, const std::vector<double> &field);

    int writeField(int bid, int zid, int sid, const std::string &fieldname, const std::vector<Vector2D> &field);

protected:

    int _fid;
};


#endif

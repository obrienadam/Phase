#ifndef PHASE_CGNS_FILE_H
#define PHASE_CGNS_FILE_H

#include <string>
#include <vector>

#include "2D/Geometry/Point2D.h"
#include "3D/Geometry/Point3D.h"

class CgnsFile
{
public:

    enum Mode
    {
        READ, WRITE, MODIFY
    };

    struct Base
    {
        std::string name;
        int id, cellDim, physDim;
    };

    struct Zone
    {
        std::string name, type;
        int id, size[9];
    };

    struct Section
    {
        std::string name, type;
        int id, start, end, nbndry, parentFlag;
        std::vector<int> cptr, cind;
    };

    struct BoCo
    {
        std::string name, type, pointListType;
        int id;
        std::vector<int> pnts;
    };

    struct Solution
    {
        std::string name, location;
        int id, dataDim, dimVals[3];
    };

    template<class T>
    struct Field
    {
        std::string name, type;
        std::array<int, 3> rmin, rmax;
        int dataDim, dimVals[3];
        std::vector<T> data;
    };

    CgnsFile()
    {}

    CgnsFile(const std::string &filename, Mode mode = READ);

    ~CgnsFile();

    int open(const std::string &filename, Mode mode);

    void close();

    int createBase(const std::string &basename, int cellDim, int physDim);

    int nBases() const;

    Base readBase(int bid) const;

    //- zones
    int nZones(int bid) const;

    Zone readZone(int bid, int zid) const;

    std::tuple<int, int> writeCoordinates(int bid, int zid, const std::vector<Point2D> &coords);

    std::tuple<int, int, int> writeCoordinates(int bid, int zid, const std::vector<Point3D> &coords);

    template<class T>
    std::vector<T> readCoords(int bid, int zid) const;

    int createStructuredZone(int bid, const std::string &zonename,
                             int nNodesI, int nNodesJ,
                             int nCellsI, int nCellsJ);

    int createStructuredZone(int bid, const std::string &zonename,
                             int nNodesI, int nNodesJ, int nNodesK,
                             int nCellsI, int nCellsJ, int nCellsK);

    int createUnstructuredZone(int bid, const std::string &zonename, int nNodes, int nCells);

    std::pair<std::vector<Label>, std::vector<Label>> readUnstructuredZone();

    int nSections(int bid, int zid) const;

    Section readSection(int bid, int zid, int sid) const;

    int writeMixedElementSection(int bid, int zid, const std::string &sectionname,
                                 int start, int end, const std::vector<int> &eptr, const std::vector<int> &eind);

    int writeBarElementSection(int bid, int zid, const std::string &sectionname,
                               int start, int end, const std::vector<int> &elements);

    int nBoCos(int bid, int zid) const;

    BoCo readBoCo(int bid, int zid, int bcid) const;

    int writeBoCo(int bid, int zid, const std::string &bcname, int start, int end);

    //- linking

    void linkNode(int bid, int zid,
                  const std::string &nodename, const std::string &filename, const std::string &nameInFile);

    void linkNode(int bid, int zid, int sid,
                  const std::string &nodename, const std::string &filename, const std::string &nameInFile);

    //- Solutions
    int nSolutions(int bid, int zid) const;

    Solution readSolution(int bid, int zid, int sid) const;

    int writeSolution(int bid, int zid, const std::string &solnname);

    int writeNodeSolution(int bid, int zid, const std::string &solname);

    //- Fields
    template<class T>
    Field<T> readField(int bid, int zid, int sid, int rmin, int rmax, const std::string& fieldname);

    template<class T>
    int writeField(int bid, int zid, int sid, const std::string &fieldname, const std::vector<T> &field);

    //- Descriptors

    int nDescriptorNodes(int bid) const;

    void writeDescriptorNode(int bid, const std::string &name, const std::string &text);

    int writeDescriptorNode(int bid, int zid, int sid, const std::string &name, const std::string &text);

protected:

    int _fid;
};


#endif

#ifndef PHASE_CGNS_FILE_H
#define PHASE_CGNS_FILE_H

#include <string>

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
    int createStructuredZone(int bid, const std::string &zonename,
                             int nNodesI, int nNodesJ,
                             int nCellsI, int nCellsJ);

    int createStructuredZone(int bid, const std::string &zonename,
                             int nNodesI, int nNodesJ, int nNodesK,
                             int nCellsI, int nCellsJ, int nCellsK);

protected:

    int _fid;
};


#endif

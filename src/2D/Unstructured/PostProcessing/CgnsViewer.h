#ifndef PHASE_CGNS_VIEWER_H
#define PHASE_CGNS_VIEWER_H

#include <vector>

#include "System/CgnsFile.h"

#include "Viewer.h"

class CgnsViewer : public Viewer
{
public:

    CgnsViewer(const Input& input, const Solver& solver);

    void write(Scalar time);

protected:

    int  createBase(int fid, const std::string& name = "Case");

    int  createZone(int fid, int bid, const FiniteVolumeGrid2D& grid, const std::string& name = "Cells");

    void writeCoords(int fid, int bid, int zid, const FiniteVolumeGrid2D& grid);

    int  writeConnectivity(int fid, int bid, int zid, const FiniteVolumeGrid2D& grid);

    void writeBoundaryConnectivity(int fid, int bid, int zid, const FiniteVolumeGrid2D& grid);

    std::string path_, gridfile_, casename_;
};

#endif

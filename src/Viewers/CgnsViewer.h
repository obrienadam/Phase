#ifndef CGNS_VIEWER_H
#define CGNS_VIEWER_H

#include <vector>

#include "Viewer.h"

class CgnsViewer : public Viewer
{
public:

    CgnsViewer(const Input& input, const Communicator& comm, const Solver& solver);

    void write(Scalar solutionTime, const Communicator& comm);

protected:

    int  createBase(int fid, const std::string& name = "Case");
    int  createZone(int fid, int bid, const FiniteVolumeGrid2D& grid, const std::string& name = "Cells");

    void writeCoords(int fid, int bid, int zid, const FiniteVolumeGrid2D& grid);
    int  writeConnectivity(int fid, int bid, int zid, const FiniteVolumeGrid2D& grid);
    void writeBoundaryConnectivity(int fid, int bid, int zid, const FiniteVolumeGrid2D& grid);
    void writeImmersedBoundaries(int fid, const Solver& solver);

    void linkGrid(int fid, int bid, int zid, const Communicator &comm);

    std::string gridfile_;
};

#endif

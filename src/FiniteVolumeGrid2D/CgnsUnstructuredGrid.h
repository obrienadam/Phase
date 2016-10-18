#ifndef CGNS_UNSTRUCTURED_GRID_H
#define CGNS_UNSTRUCTURED_GRID_H

#include "FiniteVolumeGrid2D.h"
#include "Input.h"

class CgnsUnstructuredGrid : public FiniteVolumeGrid2D
{
public:

    CgnsUnstructuredGrid();
    CgnsUnstructuredGrid(const Input& input);

    void openNewMesh(const std::string& filename, const std::string& baseName);
    void closeMesh();
    int addZone(const std::string& zoneName, int nNodes, int nCells);
    void addZoneNodes(int zoneId, const std::vector<Point2D>& nodes);

    bool fileIsOpen() const { return fileIsOpen_; }

private:

    void readNodes(int fileId, int baseId, int zoneId, int nNodes, Scalar convertToMeters);
    void readElements(int fileId, int baseId, int zoneId);
    void readBoundaries(int fileId, int baseId, int zoneId);

    int fileId_, baseId_;
    bool fileIsOpen_;
};

#endif

#ifndef PHASE_CGNS_UNSTRUCTURED_GRID_H
#define PHASE_CGNS_UNSTRUCTURED_GRID_H

#include "FiniteVolumeGrid2D.h"

class CgnsUnstructuredGrid : public FiniteVolumeGrid2D
{
public:

    CgnsUnstructuredGrid();

    CgnsUnstructuredGrid(const Input &input);

    void load(const std::string& filename);

    void readPartitionData(const std::string& filename);

private:

    void readNodes(int fileId, int baseId, int zoneId, int nNodes, Scalar convertToMeters, const Point2D& origin);

    void readElements(int fileId, int baseId, int zoneId);

    void readBoundaries(int fileId, int baseId, int zoneId);
};

#endif

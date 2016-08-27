#ifndef CGNS_UNSTRUCTURED_GRID_H
#define CGNS_UNSTRUCTURED_GRID_H

#include "FiniteVolumeGrid2D.h"
#include "Input.h"

class CgnsUnstructuredGrid : public FiniteVolumeGrid2D
{
public:
    CgnsUnstructuredGrid(const Input& input);

private:

    void readNodes(int fileId, int baseId, int zoneId, int nNodes, Scalar convertToMeters);
    void readElements(int fileId, int baseId, int zoneId);
    void readBoundaries(int fileId, int baseId, int zoneId);

};

#endif

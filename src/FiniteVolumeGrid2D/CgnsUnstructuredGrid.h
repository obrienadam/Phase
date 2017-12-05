#ifndef CGNS_UNSTRUCTURED_GRID_H
#define CGNS_UNSTRUCTURED_GRID_H

#include<cgnslib.h>

#include "FiniteVolumeGrid2D.h"
#include "Input.h"

class CgnsUnstructuredGrid : public FiniteVolumeGrid2D
{
public:

    CgnsUnstructuredGrid();

    CgnsUnstructuredGrid(const Input &input);

    void loadPartitionedGrid(std::shared_ptr<Communicator> comm);

private:

    void readNodes(int fileId, int baseId, int zoneId, int nNodes, Scalar convertToMeters, const Point2D& origin);

    void readElements(int fileId, int baseId, int zoneId);

    void readBoundaries(int fileId, int baseId, int zoneId);
};

#endif

#ifndef CGNS_UNSTRUCTURED_GRID_H
#define CGNS_UNSTRUCTURED_GRID_H

#include<cgnslib.h>

#include "FiniteVolumeGrid2D.h"
#include "Input.h"

class CgnsUnstructuredGrid : public FiniteVolumeGrid2D
{
public:

    CgnsUnstructuredGrid();
    CgnsUnstructuredGrid(const Input& input);

    //- Methods for manipulating cgns files, but not for initializing a usable mesh
    void newFile(const std::string& filename, const std::string& baseName);
    void openFile(const std::string& filename);
    void closeFile();

    int addZone(const std::string& zoneName, int nNodes, int nCells);
    void addNodes(int zoneId, const std::vector<Point2D>& nodes);
    int addTriCells(int zoneId, const std::vector<cgsize_t>& cells);
    int addMixedCells(int zoneId, int nCells, const std::vector<cgsize_t>& cells);

    int addBc(int zoneId, const std::string& name, const std::vector<cgsize_t>& faces);
    int connectZones(int zoneId, const std::vector<cgsize_t>& faces, int donorZoneId, const std::vector<cgsize_t> &donorCells);

    bool fileIsOpen() const { return fileIsOpen_; }

private:

    void readNodes(int fileId, int baseId, int zoneId, int nNodes, Scalar convertToMeters);
    void readElements(int fileId, int baseId, int zoneId);
    void readBoundaries(int fileId, int baseId, int zoneId);

    int fileId_, baseId_;
    bool fileIsOpen_;
};

#endif

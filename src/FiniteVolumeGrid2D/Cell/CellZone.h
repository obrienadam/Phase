#ifndef CELL_ZONE_H
#define CELL_ZONE_H

#include "CellGroup.h"

class CellZone : public CellGroup
{
public:

    CellZone(const std::string &name, FiniteVolumeGrid2D &grid) : CellGroup(name, grid)
    {}

    ~CellZone();

    void push_back(Cell &cell); // does not insert if the cell is found in the registry
    void moveToGroup(Cell &cell); // inserts or moves the cell to this group
    void moveToGroup(const std::vector<Ref<Cell>> &cells);

    void remove(const Cell &cell);

    void merge(CellZone &other);

    void clear();

private:

    static std::map<Label, Ref<CellZone> > registry_;
};

#endif

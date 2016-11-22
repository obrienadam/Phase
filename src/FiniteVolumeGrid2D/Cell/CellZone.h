#ifndef CELL_ZONE_H
#define CELL_ZONE_H

#include "CellGroup.h"

class CellZone : public CellGroup
{
public:

    CellZone(const std::string& name = "N/A") : CellGroup(name) {}
    ~CellZone();

    virtual void push_back(const Cell &cell); // does not insert if the cell is found in the registry
    virtual void moveToGroup(const Cell &cell); // inserts or moves the cell to this group
    virtual void moveToGroup(const std::vector<Ref<const Cell>>& cells);
    virtual void moveAllCellsToThisGroup();

    virtual void remove(const Cell &cell);
    virtual void clear();

private:

    static std::map< const Cell*, Ref<CellZone> > registry_;
};

#endif

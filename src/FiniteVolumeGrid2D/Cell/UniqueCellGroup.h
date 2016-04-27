#ifndef UNIQUE_CELL_GROUP_H
#define UNIQUE_CELL_GROUP_H

#include "CellGroup.h"

class UniqueCellGroup : public CellGroup
{
public:

    UniqueCellGroup(const std::string& name = "N/A") : CellGroup(name) {}

    virtual void push_back(const Cell &cell); // does not insert if the cell is found in the registry
    virtual void moveToGroup(const Cell &cell); // inserts or moves the cell to this group
    virtual void moveAllCellsToThisGroup();

    virtual void remove(const Cell &cell);
    virtual void clear();

private:

    static std::map< const Cell*, Ref<UniqueCellGroup> > registry_;
};

#endif

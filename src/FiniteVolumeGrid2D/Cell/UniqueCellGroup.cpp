#include "UniqueCellGroup.h"

// Initialize the registry
std::map< const Cell*, Ref<UniqueCellGroup> > UniqueCellGroup::registry_;

//- public static methods

void UniqueCellGroup::push_back(const Cell &cell)
{
    if(registry_.insert(std::make_pair(&cell, std::ref(*this))).second)
        CellGroup::push_back(cell);
}

void UniqueCellGroup::moveToGroup(const Cell &cell)
{
    auto insertion = registry_.insert(std::make_pair(&cell, std::ref(*this)));

    if(insertion.second)
        CellGroup::push_back(cell);
    else
    {
        UniqueCellGroup &group = (insertion.first)->second; // get the cell group that currently contains the cell
        group.CellGroup::remove(cell);
        CellGroup::push_back(cell);
        (insertion.first)->second = std::ref(*this);
    }
}

void UniqueCellGroup::moveAllCellsToThisGroup()
{
    for(auto &entry: registry_) // careful with this iterator
    {
        UniqueCellGroup &group = entry.second;

        if(&group == this)
            continue;

        group.::CellGroup::remove(*(entry.first));
        CellGroup::push_back(*(entry.first));
        entry.second = std::ref(*this);
    }
}

void UniqueCellGroup::remove(const Cell &cell)
{
    registry_.erase(&cell);
    CellGroup::remove(cell);
}

void UniqueCellGroup::clear()
{
    for(const Cell &cell: cells_)
        registry_.erase(&cell);

    CellGroup::clear();
}

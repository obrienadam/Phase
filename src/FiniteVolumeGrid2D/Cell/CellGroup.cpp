#include <algorithm>

#include "CellGroup.h"

void CellGroup::clear()
{
    cellSet_.clear();
    cells_.clear();
}

void CellGroup::push_back(const Cell &cell)
{
    if(cellSet_.insert(std::make_pair(&cell, cells_.size())).second)
        cells_.push_back(std::cref(cell));
}

void CellGroup::remove(const Cell &cell)
{
    auto entry = cellSet_.find(&cell);

    if(entry != cellSet_.end())
    {
        cells_.erase(cells_.begin() + entry->second);
        cellSet_.erase(entry);
    }
}

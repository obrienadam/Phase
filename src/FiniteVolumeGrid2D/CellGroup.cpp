#include "CellGroup.h"

void CellGroup::clear()
{
    directory_.clear();
    cells_.clear();
}

bool CellGroup::isInGroup(const Cell &cell)
{
    return directory_.find(&cell) == directory_.end() ? false : true;
}

void CellGroup::push_back(const Cell &cell)
{
    if(isInGroup(cell))
        return;

    cells_.push_back(Ref<const Cell>(cell));
    directory_[&cell] = cells_.size() - 1;
}

void CellGroup::add(const std::vector<Ref<const Cell> > &cells)
{
    cells_.reserve(cells.size());
    for(const Cell& cell: cells)
        push_back(cell);
}

void CellGroup::remove(const Cell &cell)
{
    if(!isInGroup(cell))
        return;

    cells_.erase(cells_.begin() + directory_[&cell]);
    directory_.erase(&cell);
}

void CellGroup::remove(const std::vector<Ref<const Cell> > &cells)
{
    for(const Cell& cell: cells)
        remove(cell);
}

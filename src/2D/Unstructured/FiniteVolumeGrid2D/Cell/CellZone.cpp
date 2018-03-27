#include "CellZone.h"

//- public static methods

CellZone::CellZone(const std::string &name,
                   const std::shared_ptr<ZoneRegistry> &registry)
        :
        Group(name)
{
    registry_ = registry;
}

CellZone::~CellZone()
{
    clear();
}

void CellZone::add(const Cell &cell)
{
    auto insertion = registry_->insert(std::make_pair(cell.id(), std::ref(*this)));

    if(!insertion.second)
    {
        insertion.first->second.get().CellGroup::remove(cell);
        insertion.first->second = std::ref(*this);
    }

    CellGroup::add(cell);
}

void CellZone::add(const CellGroup &cells)
{
    //- Must construct a temporary container since moving cells will modify the original
    //  container
    items_.reserve(size() + cells.size());
    itemSet_.reserve(size() + cells.size());

    for (const Cell &cell: std::vector<Ref<const Cell>>(cells.items()))
    {
        auto insert = registry_->insert(std::make_pair(cell.id(), std::ref(*this)));

        if(!insert.second && this != &insert.first->second.get())
        {
            insert.first->second.get().remove(cells);
            registry_->insert(std::make_pair(cell.id(), std::ref(*this)));
        }

        CellGroup::add(cell);
    }
}

void CellZone::remove(const Cell &cell)
{
    registry_->erase(cell.id());
    CellGroup::remove(cell);
}

void CellZone::remove(const CellGroup& cells)
{
    for(const Cell& cell: cells)
        if(isInGroup(cell))
            registry_->erase(cell.id());

    CellGroup::remove(cells);
}

void CellZone::clear()
{
    for (const Cell &cell: items_)
        registry_->erase(cell.id());

    CellGroup::clear();
}
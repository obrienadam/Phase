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

    if (insertion.second)
        CellGroup::add(cell);
    else
    {
        insertion.first->second.get().CellGroup::remove(cell);
        insertion.first->second = std::ref(*this);
        CellGroup::add(cell);
    }
}

void CellZone::add(const CellGroup &cells)
{
    for(const Cell& cell: cells)
        add(cell);
}

void CellZone::add(const CellZone& cells)
{
    //- Must construct a temporary container since moving cells will modify the original
    //  container

    for(const Cell &cell: cells.items())
        add(cell);
}

void CellZone::remove(const Cell &cell)
{
    if(isInGroup(cell))
    {
        registry_->erase(cell.id());
        CellGroup::remove(cell);
    }
}

void CellZone::clear()
{
    for (const Cell &cell: items_)
        registry_->erase(cell.id());

    CellGroup::clear();
}
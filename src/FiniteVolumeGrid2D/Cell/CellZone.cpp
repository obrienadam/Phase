#include "CellZone.h"

//- public static methods

CellZone::CellZone(const std::string &name,
                   const std::shared_ptr<ZoneRegistry> &registry)
    :
      Group(name)
{
    if(registry)
        registry_ = registry;
    else
        registry_ = std::shared_ptr<ZoneRegistry>(new ZoneRegistry());
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
        insertion.first->second.get().remove(cell); //- This will remove the cell from the registry
        add(cell); //- Recursive call
    }
}

void CellZone::add(const CellGroup &cells)
{
    //- Must construct a temporary container since moving cells will modify the original
    //  container
    for(const Cell& cell: cells.items())
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

void CellZone::setRegistry(std::shared_ptr<CellZone::ZoneRegistry> &registry)
{
    auto items = this->items();
    clear();
    registry_ = registry;
    CellGroup::add(items.begin(), items.end()); //- Since add(cell) is overloaded, this should work
}

#include "CellZone.h"

// Initialize the registry
std::map<Label, Ref<CellZone> > CellZone::registry_;

//- public static methods

CellZone::~CellZone()
{
    clear();
}

void CellZone::push_back(const Cell &cell)
{
    if (registry_.insert(std::make_pair(cell.id(), std::ref(*this))).second)
        CellGroup::push_back(cell);
}

void CellZone::moveToGroup(const Cell &cell)
{
    auto insertion = registry_.insert(std::make_pair(cell.id(), std::ref(*this)));

    if (insertion.second)
        CellGroup::push_back(cell);
    else
    {
        CellZone &group = (insertion.first)->second; // get the cell group that currently contains the cell
        group.CellGroup::remove(cell);
        CellGroup::push_back(cell);
        (insertion.first)->second = std::ref(*this);
    }
}

void CellZone::moveToGroup(const std::vector<Ref<Cell> > &cells)
{
    for (Cell &cell: cells)
        moveToGroup(cell);
}

void CellZone::remove(const Cell &cell)
{
    registry_.erase(cell.id());
    CellGroup::remove(cell);
}

void CellZone::merge(CellZone &other)
{
    for (const Cell &cell: other.cells())
        moveToGroup(cell);
}

void CellZone::clear()
{
    for (const Cell &cell: cells_)
        registry_.erase(cell.id());

    CellGroup::clear();
}

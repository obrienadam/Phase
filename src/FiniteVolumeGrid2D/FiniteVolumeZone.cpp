#include "FiniteVolumeZone.h"

FiniteVolumeZone::FiniteVolumeZone(const std::string &name, const std::shared_ptr<std::unordered_map<Label, Ref<CellZone> > > &cellZoneRegistry)
    :
      name_(name),
      cells_(name + "Cells", cellZoneRegistry),
      faces_(name + "Faces")
{

}

void FiniteVolumeZone::add(const Cell &cell)
{
    cells_.add(cell);

    for(const InteriorLink& nb: cell.neighbours())
        faces_.add(nb.face());

    for(const BoundaryLink& bd: cell.boundaries())
        faces_.add(bd.face());
}

void FiniteVolumeZone::add(const CellGroup &cells)
{
    for(const Cell& cell: cells)
        add(cell);
}

void FiniteVolumeZone::remove(const Cell &cell)
{
    cells_.remove(cell);

    for(const InteriorLink& nb: cell.neighbours())
        if(!cells_.isInGroup(nb.cell())) // check if a neighbour is also in the group
            faces_.remove(nb.face());

    for(const BoundaryLink& bd: cell.boundaries()) // All boundary faces get removed
        faces_.remove(bd.face());
}

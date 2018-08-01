#include "ImmersedBoundaryObject.h"
#include "FiniteVolumeGrid2D/Cell/CellSet.h"

template<class const_iterator>
std::vector<Ref<const Cell>> ImmersedBoundaryObject::internalPerimeterCells(const_iterator first, const_iterator last, bool includeDiagonals) const
{
    CellSet pCells;
    pCells.reserve(last - first);

    for(auto it = first; it != last; ++it)
    {
        if(!isInIb(*it))
            continue;

        if(includeDiagonals)
        {
            for(const CellLink &nb: it->get().cellLinks())
                if(!isInIb(nb.cell()))
                {
                    pCells.add(*it);
                    break;
                }
        }
        else
        {
            for(const CellLink &nb: it->get().neighbours())
                if(!isInIb(nb.cell()))
                {
                    pCells.add(*it);
                    break;
                }
        }
    }

    return pCells.items();
}

template<class const_iterator>
std::vector<Ref<const Cell>> ImmersedBoundaryObject::outerPerimeterCells(const_iterator first, const_iterator last, bool includeDiagonals) const
{
    CellSet pCells;
    pCells.reserve(last - first);

    for(auto it = first; it != last; ++it)
    {
        if(isInIb(*it))
            continue;

        if(includeDiagonals)
        {
            for(const CellLink &nb: it->get().cellLinks())
                if(isInIb(nb.cell()))
                {
                    pCells.add(*it);
                    break;
                }
        }
        else
        {
            for(const CellLink &nb: it->get().neighbours())
                if(isInIb(nb.cell()))
                {
                    pCells.add(*it);
                    break;
                }
        }
    }

    return pCells.items();
}

template<class T>
void ImmersedBoundaryObject::addBoundaryCondition(const std::string &name, const std::string &bcType, T bcRefValue)
{
    if(bcType == "fixed")
        addBoundaryCondition(name, ImmersedBoundaryObject::FIXED, bcRefValue);
    else if(bcType == "normal_gradient")
        addBoundaryCondition(name, ImmersedBoundaryObject::NORMAL_GRADIENT, bcRefValue);
    else if(bcType == "velocity")
        addBoundaryCondition(name, ImmersedBoundaryObject::VELOCITY, bcRefValue);
    else
        throw Exception("ImmersedBoundaryObject", "addBoundaryCondition", "invalid boundary condition type \"" + bcType + "\"");
}

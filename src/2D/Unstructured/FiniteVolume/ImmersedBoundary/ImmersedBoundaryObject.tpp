#include "ImmersedBoundaryObject.h"
#include "FiniteVolumeGrid2D/Cell/CellSet.h"

template<class const_iterator>
std::vector<Ref<const Cell>> ImmersedBoundaryObject::internalPerimeterCells(const_iterator first, const_iterator last, bool includeDiagonals) const
{
    CellSet pCells;
    pCells.reserve(last - first);

    for(auto it = first; it != last; ++it)
    {
        bool cellInIb = isInIb(*it);

        if(includeDiagonals)
            for(const CellLink &nb: it->get().cellLinks())
            {
                bool nbInIb = isInIb(nb.cell());

                if(cellInIb && !nbInIb)
                    pCells.add(*it);
                else if(!cellInIb && nbInIb)
                    pCells.add(nb.cell());
            }
        else
            for(const CellLink &nb: it->get().neighbours())
            {
                bool nbInIb = isInIb(nb.cell());

                if(cellInIb && !nbInIb)
                    pCells.add(*it);
                else if(!cellInIb && nbInIb)
                    pCells.add(nb.cell());
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
        bool cellInIb = isInIb(*it);

        if(includeDiagonals)
            for(const CellLink &nb: it->get().cellLinks())
            {
                bool nbInIb = isInIb(nb.cell());

                if(cellInIb && !nbInIb)
                    pCells.add(nb.cell());
                else if(!cellInIb && nbInIb)
                    pCells.add(*it);
            }
        else
            for(const CellLink &nb: it->get().neighbours())
            {
                bool nbInIb = isInIb(nb.cell());

                if(cellInIb && !nbInIb)
                    pCells.add(nb.cell());
                else if(!cellInIb && nbInIb)
                    pCells.add(*it);
            }
    }

    return pCells.items();
}

template<class const_iterator>
void ImmersedBoundaryObject::addIbCells(const_iterator first, const_iterator last)
{

}

template<class const_iterator>
void ImmersedBoundaryObject::addSolidCells(const_iterator first, const_iterator last)
{

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

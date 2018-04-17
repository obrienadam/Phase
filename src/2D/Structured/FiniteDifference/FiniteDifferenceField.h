#ifndef PHASE_FINITE_DIFFERENCE_FIELD_H
#define PHASE_FINITE_DIFFERENCE_FIELD_H

#include <memory.h>
#include <vector>
#include <unordered_map>
#include <string>

#include "Structured/StructuredGrid2D/StructuredGrid2D.h"

template<class T>
class FiniteDifferenceField
{
public:

    enum BcType
    {
        DIRICHLET, NEUMANN
    };

    FiniteDifferenceField(const std::shared_ptr<const StructuredGrid2D> &grid);

    T &operator()(size_t i, size_t j)
    { return _nodes[j * _grid->nNodesI() + i]; }

    T operator()(size_t i, size_t j) const
    { return _nodes[j * _grid->nNodesI() + i]; }

protected:

    std::shared_ptr<const StructuredGrid2D> _grid;

    std::vector<T> _nodes;

    std::unordered_map<StructuredGrid2D::Boundary, std::pair<BcType, T>> _bcs;
};

#include "FiniteDifferenceField.tpp"

#endif //PHASE_FIELD_H

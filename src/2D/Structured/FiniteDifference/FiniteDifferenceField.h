#ifndef PHASE_FINITE_DIFFERENCE_FIELD_H
#define PHASE_FINITE_DIFFERENCE_FIELD_H

#include <memory.h>
#include <vector>
#include <unordered_map>
#include <string>

#include "Structured/StructuredGrid2D/StructuredGrid2D.h"
#include "IndexMap.h"

template<class T>
class FiniteDifferenceField
{
public:

    enum BcType
    {
        DIRICHLET, NEUMANN
    };

    FiniteDifferenceField(const std::shared_ptr<const StructuredGrid2D> &grid, const std::string &name = "");

    const std::shared_ptr<const StructuredGrid2D> &grid() const
    { return _grid; }

    constexpr const std::string &name() const
    { return _name; }

    const std::shared_ptr<IndexMap> &idxMap() const
    { return _idxMap; }

    T &operator()(size_t i, size_t j)
    { return _nodes[j * _grid->nNodesI() + i]; }

    T operator()(size_t i, size_t j) const
    { return _nodes[j * _grid->nNodesI() + i]; }

protected:

    std::string _name;

    std::shared_ptr<const StructuredGrid2D> _grid;

    std::vector<T> _nodes;

    std::unordered_map<StructuredGrid2D::Boundary, std::pair<BcType, T>> _bcs;

    std::shared_ptr<IndexMap> _idxMap;
};

#include "FiniteDifferenceField.tpp"

#endif //PHASE_FIELD_H

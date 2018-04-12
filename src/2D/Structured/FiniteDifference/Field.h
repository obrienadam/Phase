#ifndef PHASE_FIELD_H
#define PHASE_FIELD_H

#include <memory.h>
#include <vector>

#include "StructuredGrid2D/StructuredGrid2D.h"

template<class T>
class Field
{
public:

    Field(const std::shared_ptr<const StructuredGrid2D> &grid, const std::string &name = "N/A")
            :
            std::vector<T>::std::vector(grid->nNodes()),
            grid_(grid),
            name(name)
    {}

    T &operator()(size_t i, size_t j)
    { return std::vector<T>::operator[](grid_->id(i, j)); }

    const T &operator()(size_t i, size_t j) const
    { return std::vector<T>::operator[](grid_->id(i, j)); }

    std::string name;

protected:

    std::vector<T> nodevals_, ifacevals_, jfacevals_;

    std::shared_ptr<const StructuredGrid2D> grid_;
};


#endif //PHASE_FIELD_H

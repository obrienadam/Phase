#ifndef CELL_GROUP_H
#define CELL_GROUP_H

#include <vector>
#include <map>

#include "Types.h"
#include "Cell.h"

class CellGroup
{
public:

    CellGroup(const std::string& name = "N/A") : name_(name) {}

    const std::string& name() const { return name_; }
    const std::string& rename(const std::string& name) { return name_ = name; }

    void clear();
    void reserve(size_t size) { cells_.reserve(size); }
    size_t size() const { return cells_.size(); }

    const std::vector< Ref<const Cell> >& cells() const { return cells_; }

    bool isInGroup(const Cell &cell);
    void push_back(const Cell &cell);
    void add(const std::vector< Ref<const Cell> >& cells);

    void remove(const Cell &cell);
    void remove(const std::vector< Ref<const Cell> >& cells);

    std::vector< Ref<const Cell> >::iterator begin() { return cells_.begin(); }
    std::vector< Ref<const Cell> >::iterator end() { return cells_.end(); }
    std::vector< Ref<const Cell> >::const_iterator begin() const { return cells_.begin(); }
    std::vector< Ref<const Cell> >::const_iterator end() const { return cells_.end(); }

private:

    std::string name_;

    std::map< const Cell *, size_t> directory_;
    std::vector< Ref<const Cell> > cells_;
};

#endif

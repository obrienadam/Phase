#ifndef CELL_ZONE_H
#define CELL_ZONE_H

#include "CellGroup.h"

class CellZone : public CellGroup
{
public:

    typedef std::unordered_map<Label, Ref<CellZone>> ZoneRegistry;

    CellZone(const std::string &name = "N/A",
             const std::shared_ptr<ZoneRegistry>& registry = std::make_shared<ZoneRegistry>());

    ~CellZone();

    void add(const Cell &cell);

    void add(const CellGroup& cells);

    template <class const_iterator>
    void add(const_iterator begin, const_iterator end)
    {
        for(const_iterator itr = begin; itr != end; ++itr)
            add(*itr);
    }

    void remove(const Cell &cell);

    void clear();

    std::shared_ptr<ZoneRegistry> registry() const
    { return registry_; }

private:

    std::shared_ptr<ZoneRegistry> registry_;
};

#endif

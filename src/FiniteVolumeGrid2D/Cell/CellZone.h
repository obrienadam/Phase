#ifndef CELL_ZONE_H
#define CELL_ZONE_H

#include "CellGroup.h"

class CellZone : public CellGroup
{
public:

    typedef std::unordered_map<Label, Ref<CellZone>> ZoneRegistry;

    CellZone(const std::string &name = "N/A",
             const std::shared_ptr<ZoneRegistry>& registry = nullptr);

    ~CellZone();

    void add(const Cell &cell);

    void add(const CellGroup& cells);

    void remove(const Cell &cell);

    void clear();

    void setRegistry(std::shared_ptr<ZoneRegistry>& registry);

    std::shared_ptr<ZoneRegistry> registry() const
    { return registry_; }

private:

    std::shared_ptr<ZoneRegistry> registry_;
};

#endif

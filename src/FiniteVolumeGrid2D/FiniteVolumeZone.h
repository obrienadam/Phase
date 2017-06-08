#ifndef FINITE_VOLUME_ZONE_H
#define FINITE_VOLUME_ZONE_H

#include "CellZone.h"
#include "FaceGroup.h"

class FiniteVolumeZone
{
public:

    FiniteVolumeZone(const std::string& name,
                     const std::shared_ptr<std::unordered_map<Label, Ref<CellZone>>>& cellZoneRegistry);

    const std::string& name() const
    { return name_; }

    const CellZone& cells() const
    { return cells_; }

    const FaceGroup& faces() const
    { return faces_; }

    void add(const Cell& cell);

    void add(const CellGroup& cells);

    void remove(const Cell& cell);

private:

    CellZone cells_;

    FaceGroup faces_;

    std::string name_;
};

#endif

#ifndef EQUIDISTANT_GRID_2D
#define EQUIDISTANT_GRID_2D

#include "FiniteVolumeGrid2D.h"

class EquidistantGrid2D : public FiniteVolumeGrid2D
{
public:

    EquidistantGrid2D(int nCellsI, int nCellsJ, Scalar h);

    virtual Vector2D sf(int i, int j, Face face) const;
    virtual Vector2D rf(int i, int j, Face face) const;
    virtual Vector2D rc(int i, int j, Face face) const;

    Scalar cellVolume(int i, int j) const { return h_*h_; }

private:

    Scalar h_;
};

#endif

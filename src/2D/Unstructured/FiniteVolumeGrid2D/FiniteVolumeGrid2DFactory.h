#ifndef FINITE_VOLUME_GRID_2D_FACTORY_H
#define FINITE_VOLUME_GRID_2D_FACTORY_H

#include "FiniteVolumeGrid2D.h"

class FiniteVolumeGrid2DFactory
{
public:

    enum GridType
    {
        CGNS,
        RECTILINEAR,
        RELOAD
    };

    static std::shared_ptr<FiniteVolumeGrid2D> create(GridType type, const Input &input);

    static std::shared_ptr<FiniteVolumeGrid2D> create(std::string type, const Input &input);

    static std::shared_ptr<FiniteVolumeGrid2D> create(const Input &input);
};


#endif //PHASE_FINITEVOLUMEGRID2DFACTORY_H

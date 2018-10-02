#ifndef PHASE_CELESTE_AXISYMMETRIC_IMMERSED_BOUNDARY_H
#define PHASE_CELESTE_AXISYMMETRIC_IMMERSED_BOUNDARY_H

#include "CelesteImmersedBoundary.h"

class CelesteAxisymmetricImmersedBoundary: public CelesteImmersedBoundary
{
public:

    CelesteAxisymmetricImmersedBoundary(const Input &input,
                                        const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                        const std::shared_ptr<CellGroup> &fluidCells,
                                        const std::weak_ptr<const ImmersedBoundary> &ib);



protected:

    virtual void computeCurvature() override;

    std::shared_ptr<VectorFiniteVolumeField> kappaRZ_;
};

#endif

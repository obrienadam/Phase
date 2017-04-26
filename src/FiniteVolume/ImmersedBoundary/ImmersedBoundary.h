#ifndef IMMERSED_BOUNDARY_H
#define IMMERSED_BOUNDARY_H

#include "ImmersedBoundaryObject.h"

class Solver;

class ImmersedBoundary
{
public:

    enum
    {
        FLUID = 1, IB = 2, SOLID = 3, FRESHLY_CLEARED = 4, BUFFER=5
    };

    ImmersedBoundary(const Input &input, const Communicator &comm, Solver &solver);

    void initCellZones();

    void update(Scalar timeStep);

    Size nIbObjects() const
    { return ibObjs_.size(); }

    std::vector<Ref<const ImmersedBoundaryObject>> ibObjs() const;

    bool isIbCell(const Cell &cell) const;

    void cutFaces();

protected:

    void setCellStatus();

    Solver &solver_;
    const Communicator &comm_;
    FiniteVolumeField<int> &cellStatus_;
    std::vector<std::shared_ptr<ImmersedBoundaryObject>> ibObjs_;

};

#endif

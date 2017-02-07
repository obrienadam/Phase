#ifndef IMMERSED_BOUNDARY_H
#define IMMERSED_BOUNDARY_H

#include "ImmersedBoundaryObject.h"

class Solver;

class ImmersedBoundary
{
public:

    enum {FLUID = 1, IB = 2, SOLID = 3};

    ImmersedBoundary(const Input& input, const Communicator& comm, Solver &solver);

    void initCellZones(const Communicator& comm);

    Equation<Scalar> eqns(ScalarFiniteVolumeField& field);
    Equation<Vector2D> eqns(VectorFiniteVolumeField& field);

    const std::vector<ImmersedBoundaryObject>& ibObjs() const { return ibObjs_; }
    bool isIbCell(const Cell& cell) const;

protected:

    void setCellStatus(const Communicator& comm);

    const Solver &solver_;
    ScalarFiniteVolumeField &cellStatus_;
    std::vector<ImmersedBoundaryObject> ibObjs_;

};

#endif

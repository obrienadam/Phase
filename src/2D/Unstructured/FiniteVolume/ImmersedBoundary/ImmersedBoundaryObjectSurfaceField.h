#ifndef PHASE_IMMERSED_BOUNDARY_OBJECT_SURFACE_FIELD_H
#define PHASE_IMMERSED_BOUNDARY_OBJECT_SURFACE_FIELD_H

#include "System/Communicator.h"

#include "ImmersedBoundaryObject.h"

template<class T>
class ImmersedBoundaryObject::SurfaceField
{
public:

    SurfaceField(const std::shared_ptr<const ImmersedBoundaryObject> &ibObj);

    Size size() const
    { return _field.size(); }

    const Point2D &pt(Label i) const
    { return _field[i].first; }

    T& operator()(Label i)
    { return _field[i].second; }

    const T& operator()(Label i) const
    { return _field[i].second; }

    void sortCounterClockwise();

    void gather(const Communicator &comm, int procNo);

private:

    std::shared_ptr<const ImmersedBoundaryObject> _ibObj;

    std::vector<std::pair<Point2D, T>> _field;
};

#include "ImmersedBoundaryObjectSurfaceField.tpp"

#endif

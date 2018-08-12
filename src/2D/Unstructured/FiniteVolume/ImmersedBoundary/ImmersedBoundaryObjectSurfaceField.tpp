#include "ImmersedBoundaryObjectSurfaceField.h"

template<class T>
ImmersedBoundaryObject::SurfaceField<T>::SurfaceField(const std::shared_ptr<const ImmersedBoundaryObject> &ibObj)
    :
      _ibObj(ibObj)
{
    _field.reserve(_ibObj->ibCells().size());

    for(const Cell &c: _ibObj->ibCells())
        _field.push_back(std::make_pair(_ibObj->nearestIntersect(c.centroid()), T()));
}

template<class T>
void ImmersedBoundaryObject::SurfaceField<T>::sortCounterClockwise()
{
    std::sort(_field.begin(), _field.end(), [this](const std::pair<Point2D, T> &lhs, const std::pair<Point2D, T> &rhs){
        return (lhs.first - _ibObj->shape().centroid()).angle() < (rhs.first - _ibObj->shape().centroid()).angle();
    });
}

template<class T>
void ImmersedBoundaryObject::SurfaceField<T>::gather(const Communicator &comm, int procNo)
{
    _field = comm.gatherv(procNo, _field);
}

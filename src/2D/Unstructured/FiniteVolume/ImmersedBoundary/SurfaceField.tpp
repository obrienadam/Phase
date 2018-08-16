#include "SurfaceField.h"

template<class T>
SurfaceField<T>::SurfaceField(const std::shared_ptr<const ImmersedBoundaryObject> &ibObj)
    :
      _ibObj(ibObj)
{
    _field.reserve(_ibObj->ibCells().size());

    for(const Cell &c: _ibObj->ibCells())
        _field.push_back(std::make_pair(_ibObj->nearestIntersect(c.centroid()), T()));
}

template<class T>
void SurfaceField<T>::sortCounterClockwise(const Point2D &pt)
{
    std::sort(_field.begin(), _field.end(), [&pt](const std::pair<Point2D, T> &lhs, const std::pair<Point2D, T> &rhs){
        return (lhs.first - pt).angle() < (rhs.first - pt).angle();
    });
}

template<class T>
Scalar SurfaceField<T>::integrate()
{
    Scalar result = 0.;

    for(const auto &pt: _field)
    {

    }

    return result;
}

template<class T>
void SurfaceField<T>::gather(const Communicator &comm, int procNo)
{
    _field = comm.gatherv(procNo, _field);
}

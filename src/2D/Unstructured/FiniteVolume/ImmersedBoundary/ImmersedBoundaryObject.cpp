#include "System/NotImplementedException.h"

#include "ImmersedBoundaryObject.h"
#include "ImmersedBoundary.h"

ImmersedBoundaryObject::ImmersedBoundaryObject(const std::string &name)
    :
      _name(name)
{
    _cells = CellGroup("Cells");
    _ibCells = CellGroup("IbCells");
    _solidCells = CellGroup("SolidCells");

    _force = Vector2D(0., 0.);
    _torque = 0.;
}

void ImmersedBoundaryObject::setIbCells(const CellGroup &ibCells)
{
    _cells.remove(_ibCells);
    _cells.add(ibCells);
    _ibCells = ibCells;
}

void ImmersedBoundaryObject::setSolidCells(const CellGroup &solidCells)
{
    _cells.remove(_solidCells);
    _cells.add(solidCells);
    _solidCells = solidCells;
}

bool ImmersedBoundaryObject::addIbCell(const Cell &cell)
{
    if(_cells.add(cell))
        return _ibCells.add(cell);
    return false;
}

bool ImmersedBoundaryObject::addSolidCell(const Cell &cell)
{
    if(_cells.add(cell))
        return _solidCells.add(cell);
    return false;
}

void ImmersedBoundaryObject::clear()
{
    _cells.clear();
    _ibCells.clear();
    _solidCells.clear();
}

void ImmersedBoundaryObject::initCircle(const Point2D &center, Scalar radius)
{
    _shape = std::unique_ptr<Circle>(new Circle(center, radius));
}


void ImmersedBoundaryObject::initBox(const Point2D &center, Scalar width, Scalar height)
{
    _shape = std::unique_ptr<Box>(new Box(
                                      Point2D(center.x - width / 2., center.y - height / 2.),
                                      Point2D(center.x + width / 2., center.y + height / 2.)
                                      ));
}

LineSegment2D ImmersedBoundaryObject::intersectionLine(const LineSegment2D &ln) const
{
    auto xc = _shape->intersections(ln);

    if (xc.empty())
    {
        Point2D pts[] = {ln.ptA(), ln.ptB()};
        xc.push_back(_shape->closest(pts, pts + 2));
    }

    return LineSegment2D(ln.ptA(), xc[0]);
}

LineSegment2D ImmersedBoundaryObject::intersectionLine(const Point2D &ptA, const Point2D &ptB) const
{
    return intersectionLine(LineSegment2D(ptA, ptB));
}

Vector2D ImmersedBoundaryObject::nearestEdgeUnitNormal(const Point2D &pt) const
{
    switch (_shape->type())
    {
    case Shape2D::CIRCLE:
        return (_shape->centroid() - pt).unitVec();
    case Shape2D::BOX:
    case Shape2D::POLYGON:
    {
        auto edge = _shape->nearestEdge(pt);
        return dot(edge.norm(), _shape->centroid() - edge.center()) > 0. ? edge.norm().unitVec()
                                                                         : -edge.norm().unitVec();
    }
    default:
        throw Exception("ImmersedBoundaryObject", "nearestEdgeNormal", "not implemented for specified shape.");
    }
}

ImmersedBoundaryObject::BoundaryConditionType ImmersedBoundaryObject::bcType(const std::string &name) const
{
    auto it = _bcTypes.find(name);
    return it == _bcTypes.end() ? NORMAL_GRADIENT: it->second;
}

template<>
Scalar ImmersedBoundaryObject::bcRefValue<Scalar>(const std::string &name) const
{
    auto it = _bcScalarRefValues.find(name);
    return it == _bcScalarRefValues.end() ? 0. : it->second;
}

template<>
Vector2D ImmersedBoundaryObject::bcRefValue<Vector2D>(const std::string &name) const
{
    auto it = _bcVectorRefValues.find(name);
    return it == _bcVectorRefValues.end() ? Vector2D(0., 0.) : it->second;
}

template<>
void ImmersedBoundaryObject::addBoundaryCondition(const std::string &name, BoundaryConditionType bcType, Scalar bcRefValue)
{
    _bcTypes[name] = bcType;
    _bcScalarRefValues[name] = bcRefValue;
}

template<>
void ImmersedBoundaryObject::addBoundaryCondition(const std::string &name, BoundaryConditionType bcType, Vector2D bcRefValue)
{
    _bcTypes[name] = bcType;
    _bcVectorRefValues[name] = bcRefValue;
}

template<>
void ImmersedBoundaryObject::addBoundaryCondition(const std::string &name, BoundaryConditionType bcType, std::string bcRefValue)
{
    try
    {
        addBoundaryCondition(name, bcType, std::stod(bcRefValue));
    }
    catch(const Exception &e)
    {
        throw;
    }
    catch (const std::exception &e)
    {
        try
        {
            addBoundaryCondition(name, bcType, Vector2D(bcRefValue));
        }
        catch(const std::exception &e)
        {
            throw Exception("ImmersedBoundaryObject", "addBoundaryCondition", "invalid reference value \"" + bcRefValue + "\".");
        }
    }
}

Vector2D ImmersedBoundaryObject::velocity(const Point2D &point) const
{
    return _motion ? _motion->velocity(point) : Vector2D(0., 0.);
}

Vector2D ImmersedBoundaryObject::acceleration(const Point2D &point) const
{
    return _motion ? _motion->acceleration(point) : Vector2D(0., 0.);
}

Scalar ImmersedBoundaryObject::theta() const
{
    return _motion ? _motion->theta() : 0.;
}

Scalar ImmersedBoundaryObject::omega() const
{
    return _motion ? _motion->omega() : 0.;
}

Scalar ImmersedBoundaryObject::alpha() const
{
    return _motion ? _motion->alpha() : 0.;
}

void ImmersedBoundaryObject::updatePosition(Scalar timeStep)
{
    if (_motion)
    {
        _motion->update(timeStep);
        _shape->move(_motion->position());
    }
}

#include <math.h>

#include <boost/algorithm/string.hpp>

#include "Solver.h"
#include "FaceInterpolation.h"
#include "EigenSparseMatrixSolver.h"

Solver::Solver(const Input &input, const Communicator &comm, FiniteVolumeGrid2D &grid)
    :
      ibObjManager_(input, comm, *this),
      comm_(comm),
      grid_(grid)
{
    //- Set simulation time options
    std::string timeDependentOpt = input.caseInput().get<std::string>("Solver.timeDependent");
    boost::to_lower(timeDependentOpt);

    timeDependent_ = timeDependentOpt == "on" ? ON : OFF;
    timeStepRelaxation_ = input.caseInput().get<Scalar>("Solver.timeStepRelaxation", 1.);
    maxTimeStep_ = input.caseInput().get<Scalar>("Solver.timeStep");

    //- Set voluem integrators
    volumeIntegrators_ = VolumeIntegrator::initVolumeIntegrators(input, *this);

    int rank = comm.rank();

    FiniteVolumeField<int>& proc = addIntegerField("proc");

    for(const Cell& cell: grid_.cells())
        proc(cell) = rank;
}

std::string Solver::info() const
{
    return "SOLVER INFO\n"
           "Time dependent: " + std::string((timeDependent_ == ON) ? "On" : "Off") + "\n";
}

FiniteVolumeField<int> &Solver::addIntegerField(const std::string &name)
{
    auto insert = integerFields_.insert(std::make_pair(name, FiniteVolumeField<int>(grid_, name)));

    if(!insert.second)
        throw Exception("Solver", "addIntegerField", "field \"" + name + "\" already exists.");

    return insert.first->second;
}

ScalarFiniteVolumeField& Solver::addScalarField(const Input& input, const std::string& name)
{
    auto insert = scalarFields_.insert(std::make_pair(name, ScalarFiniteVolumeField(input, grid_, name)));

    if(!insert.second)
        throw Exception("Solver", "addScalarField", "field \"" + name + "\" already exists.");


    return insert.first->second;
}

ScalarFiniteVolumeField& Solver::addScalarField(const std::string& name)
{
    auto insert = scalarFields_.insert(std::make_pair(name, ScalarFiniteVolumeField(grid_, name)));

    if(!insert.second)
        throw Exception("Solver", "addScalarField", "field \"" + name + "\" already exists.");

    return insert.first->second;
}

VectorFiniteVolumeField& Solver::addVectorField(const Input& input, const std::string& name)
{
    auto insert = vectorFields_.insert(std::make_pair(name, VectorFiniteVolumeField(input, grid_, name)));

    if(!insert.second)
        throw Exception("Solver", "addVectorField", "field \"" + name + "\" already exists.");

    return insert.first->second;
}

VectorFiniteVolumeField& Solver::addVectorField(const std::string& name)
{
    auto insert = vectorFields_.insert(std::make_pair(name, VectorFiniteVolumeField(grid_, name)));

    if(!insert.second)
        throw Exception("Solver", "addVectorField", "field \"" + name + "\" already exists.");

    return insert.first->second;
}

std::vector<Polygon>& Solver::addGeometries(const std::string &name)
{
    return (geometries_.insert(std::make_pair(name, std::vector<Polygon>(grid_.cells().size()))).first)->second;
}

void Solver::setInitialConditions(const Input& input)
{
    using namespace std;
    using namespace boost::property_tree;

    for(const auto& child: input.initialConditionInput().get_child("InitialConditions"))
    {
        auto scalarFieldIt = scalarFields().find(child.first);

        if(scalarFieldIt != scalarFields().end())
        {
            ScalarFiniteVolumeField &field = scalarFieldIt->second;

            for(const auto& ic: child.second)
            {
                const auto &icTree = ic.second;
                std::string type = icTree.get<string>("type");

                if(type == "circle")
                {
                    Circle circle = Circle(Vector2D(icTree.get<string>("center")), icTree.get<Scalar>("radius"));
                    setCircle(circle, icTree.get<Scalar>("value"), field);
                }
                else if(type == "circleSector")
                {
                    Circle circle = Circle(Vector2D(icTree.get<string>("center")), icTree.get<Scalar>("radius"));
                    Scalar thetaMin = icTree.get<Scalar>("thetaMin");
                    Scalar thetaMax = icTree.get<Scalar>("thetaMax");
                    Scalar innerValue = icTree.get<Scalar>("value");
                    setCircleSector(circle, thetaMin, thetaMax, innerValue, field);
                }
                else if(type == "box")
                {
                    Point2D center = Point2D(icTree.get<string>("center"));
                    Scalar w = icTree.get<Scalar>("width")/2;
                    Scalar h = icTree.get<Scalar>("height")/2;

                    std::vector<Point2D> vertices = {
                        Point2D(center.x - w, center.y - h),
                        Point2D(center.x + w, center.y - h),
                        Point2D(center.x + w, center.y + h),
                        Point2D(center.x - w, center.y + h)
                    };

                    setBox(Polygon(vertices), icTree.get<Scalar>("value"), field);
                }
                else if(type == "uniform")
                    field.fillInterior(icTree.get<Scalar>("value"));
                else if(type == "rotating")
                {
                    setRotating(icTree.get<std::string>("function"),
                                icTree.get<Scalar>("amplitude"),
                                Vector2D(icTree.get<std::string>("center")),
                                field);
                }
                else
                    throw Exception("Input", "setInitialConditions", "invalid initial condition type \"" + type + "\".");

                printf("Set initial condition \"%s\" of type %s on field \"%s\".\n", ic.first.c_str(), type.c_str(), field.name().c_str());
            }

            continue;
        }

        auto vectorFieldIt = vectorFields().find(child.first);

        if(vectorFieldIt != vectorFields().end())
        {
            VectorFiniteVolumeField &field = vectorFieldIt->second;

            for(const auto& ic: child.second)
            {
                const auto &icTree = ic.second;
                std::string type = icTree.get<string>("type");

                if(type == "circle")
                {
                    Circle circle = Circle(Vector2D(icTree.get<string>("center")), icTree.get<Scalar>("radius"));
                    setCircle(circle, Vector2D(icTree.get<string>("value")), field);
                }
                else if(type == "square")
                {
                    Point2D center = Point2D(icTree.get<string>("center"));
                    Scalar w = icTree.get<Scalar>("width")/2;
                    Scalar h = icTree.get<Scalar>("height")/2;

                    std::vector<Point2D> vertices = {
                        Point2D(center.x - w, center.y - h),
                        Point2D(center.x + w, center.y - h),
                        Point2D(center.x + w, center.y + h),
                        Point2D(center.x - w, center.y + h)
                    };

                    setBox(Polygon(vertices), Vector2D(icTree.get<string>("value")), field);
                }
                else if(type == "uniform")
                    field.fillInterior(Vector2D(icTree.get<string>("value")));
                else if(type == "rotating")
                {
                    setRotating(icTree.get<std::string>("xFunction"),
                                icTree.get<std::string>("yFunction"),
                                Vector2D(icTree.get<std::string>("amplitude")),
                                Vector2D(icTree.get<std::string>("center")),
                                field);
                }
                else
                    throw Exception("Input", "setInitialConditions", "invalid initial condition type \"" + type + "\".");

                printf("Set initial condition \"%s\" of type %s on field \"%s\".\n", ic.first.c_str(), type.c_str(), field.name().c_str());
            }
        }
    }
}

//- Protected methods

void Solver::setCircle(const Circle &circle, Scalar innerValue, ScalarFiniteVolumeField &field)
{
    const Polygon pgn = circle.polygonize(1000);

    for(const Cell& cell: field.grid.cells())
    {
        const Polygon xc = intersectionPolygon(pgn, cell.shape());

        if(xc.area() == 0.)
            continue;

        field(cell) = innerValue*xc.area()/cell.shape().area();
    }

    interpolateFaces(fv::INVERSE_VOLUME, field);
}

void Solver::setCircle(const Circle &circle, const Vector2D &innerValue, VectorFiniteVolumeField &field)
{
    const Polygon pgn = circle.polygonize(400);

    for(const Cell& cell: field.grid.cells())
    {
        const Polygon xc = intersectionPolygon(pgn, cell.shape());

        if(xc.area() == 0.)
            continue;

        field(cell) = innerValue*xc.area()/cell.shape().area();
    }

    interpolateFaces(fv::INVERSE_VOLUME, field);
}

void Solver::setCircleSector(const Circle &circle, Scalar thetaMin, Scalar thetaMax, Scalar innerValue, ScalarFiniteVolumeField &field)
{
    thetaMin *= M_PI/180;
    thetaMax *= M_PI/180;

    while(thetaMin > 2*M_PI)
        thetaMin -= 2*M_PI;

    while(thetaMin < 0.)
        thetaMin += 2*M_PI;

    while(thetaMax > 2*M_PI)
        thetaMax -= 2*M_PI;

    while(thetaMax < 0.)
        thetaMax += 2*M_PI;

    std::vector<Point2D> vtx;
    const int nPts = 400;
    const Scalar dTheta = (thetaMax - thetaMin)/(nPts - 1);

    for(int i = 0; i < nPts; ++i)
    {
        const Scalar theta = thetaMin + i*dTheta;
        const Vector2D rVec(circle.radius()*cos(theta), circle.radius()*sin(theta));

        vtx.push_back(circle.centroid() + rVec);
    }

    const Polygon pgn(vtx);

    for(const Cell& cell: field.grid.cells())
    {
        const Polygon xc = intersectionPolygon(pgn, cell.shape());

        if(xc.area() == 0.)
            continue;

        field(cell) = innerValue*xc.area()/cell.shape().area();
    }

    interpolateFaces(fv::INVERSE_VOLUME, field);
}

void Solver::setBox(const Polygon& box, Scalar innerValue, ScalarFiniteVolumeField& field)
{
    for(const Cell& cell: field.grid.cells())
    {
        const Polygon xc = intersectionPolygon(box, cell.shape());

        if(xc.area() == 0.)
            continue;

        field(cell) = innerValue*xc.area()/cell.shape().area();
    }

    fv::interpolateFaces(fv::INVERSE_VOLUME, field);
}

void Solver::setBox(const Polygon& box, const Vector2D& innerValue, VectorFiniteVolumeField& field)
{
    for(const Cell& cell: field.grid.cells())
    {
        const Polygon xc = intersectionPolygon(box, cell.shape());

        if(xc.area() == 0.)
            continue;

        field(cell) = innerValue*xc.area()/cell.shape().area();
    }

    fv::interpolateFaces(fv::INVERSE_VOLUME, field);
}

void Solver::setRotating(const std::string &function, Scalar amplitude, const Vector2D &center, ScalarFiniteVolumeField &field)
{
    std::function<Scalar(Scalar)> func;

    if(function == "sin")
        func = [](Scalar x){ return sin(x); };
    else if(function == "cos")
        func = [](Scalar x){ return cos(x); };
    else
        throw Exception("Input", "setRotating", "invalid rotation function.");

    for(const Cell& cell: field.grid.cells())
    {
        Vector2D rVec = cell.centroid() - center;
        Scalar theta = atan2(rVec.y, rVec.x);

        field[cell.id()] = amplitude*func(theta);
    }

    for(const Face& face: field.grid.interiorFaces())
    {
        Vector2D rVec = face.centroid() - center;
        Scalar theta = atan2(rVec.y, rVec.x);

        field.faces()[face.id()] = amplitude*func(theta);
    }
}

void Solver::setRotating(const std::string &xFunction, const std::string &yFunction, const Vector2D &amplitude, const Vector2D &center, VectorFiniteVolumeField &field)
{
    std::function<Scalar(Scalar)> xFunc, yFunc;

    if(xFunction == "sin")
        xFunc = [](Scalar x) { return sin(x); };
    else if(xFunction == "cos")
        xFunc = [](Scalar x) { return cos(x); };
    else
        throw Exception("Input", "setRotating", "invalid x rotation function.");

    if(yFunction == "sin")
        yFunc = [](Scalar x) { return sin(x); };
    else if(yFunction == "cos")
        yFunc = [](Scalar x) { return cos(x); };
    else
        throw Exception("Input", "setRotating", "invalid y rotation function.");

    for(const Cell& cell: field.grid.cells())
    {
        Vector2D rVec = cell.centroid() - center;
        Scalar theta = atan2(rVec.y, rVec.x);

        field[cell.id()].x = amplitude.x*xFunc(theta);
        field[cell.id()].y = amplitude.y*yFunc(theta);
    }

    for(const Face& face: field.grid.interiorFaces())
    {
        Vector2D rVec = face.centroid() - center;
        Scalar theta = atan2(rVec.y, rVec.x);

        field.faces()[face.id()].x = amplitude.x*xFunc(theta);
        field.faces()[face.id()].y = amplitude.y*yFunc(theta);
    }
}

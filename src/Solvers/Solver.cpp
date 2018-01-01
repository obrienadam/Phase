#include <math.h>
#include <regex>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <cgnslib.h>

#include "Solver.h"
#include "FaceInterpolation.h"
#include "EigenSparseMatrixSolver.h"

Solver::Solver(const Input &input, std::shared_ptr<FiniteVolumeGrid2D> &grid)
        :
        ib_(input, *this),
        grid_(grid)
{
    //- Set simulation time options
    maxTimeStep_ = input.caseInput().get<Scalar>("Solver.timeStep");
}

void Solver::printf(const char *format, ...) const
{
    va_list argsPtr;
    va_start(argsPtr, format);

    if (grid_->comm().isMainProc())
        vfprintf(stdout, format, argsPtr);

    va_end(argsPtr);
}

Scalar Solver::getStartTime(const Input &input) const
{
    if (input.initialConditionInput().get<std::string>("InitialConditions.type", "") == "restart")
    {
        using namespace std;
        using namespace boost::filesystem;

        std::regex re("[0-9]+\\.[0-9]+");
        Scalar maxTime = 0.;

        for (directory_iterator end, dir("./solution"); dir != end; ++dir)
            if (regex_match(dir->path().filename().string(), re))
            {
                std::smatch match;
                std::regex_search(dir->path().filename().string(), match, re);
                Scalar time = std::stod(match.str());

                maxTime = std::max(time, maxTime);
            }

        return maxTime;
    }

    return 0.;
}

FiniteVolumeField<int> &Solver::addIntegerField(const std::string &name)
{
    auto insert = integerFields_.insert(std::make_pair(name, std::make_shared<FiniteVolumeField < int>>
    (grid_, name)));

    if (!insert.second)
        throw Exception("Solver", "addIntegerField", "field \"" + name + "\" already exists.");

    return *insert.first->second;
}

ScalarFiniteVolumeField &Solver::addScalarField(const Input &input, const std::string &name)
{
    auto insert = scalarFields_.insert(
            std::make_pair(name, std::make_shared<ScalarFiniteVolumeField>(input, grid_, name)));

    if (!insert.second)
        throw Exception("Solver", "addScalarField", "field \"" + name + "\" already exists.");


    return *insert.first->second;
}

ScalarFiniteVolumeField &Solver::addScalarField(const std::string &name)
{
    auto insert = scalarFields_.insert(std::make_pair(name, std::make_shared<ScalarFiniteVolumeField>(grid_, name)));

    if (!insert.second)
        throw Exception("Solver", "addScalarField", "field \"" + name + "\" already exists.");

    return *insert.first->second;
}

VectorFiniteVolumeField &Solver::addVectorField(const Input &input, const std::string &name)
{
    auto insert = vectorFields_.insert(
            std::make_pair(name, std::make_shared<VectorFiniteVolumeField>(input, grid_, name)));

    if (!insert.second)
        throw Exception("Solver", "addVectorField", "field \"" + name + "\" already exists.");

    return *insert.first->second;
}

VectorFiniteVolumeField &Solver::addVectorField(const std::string &name)
{
    auto insert = vectorFields_.insert(std::make_pair(name, std::make_shared<VectorFiniteVolumeField>(grid_, name)));

    if (!insert.second)
        throw Exception("Solver", "addVectorField", "field \"" + name + "\" already exists.");

    return *insert.first->second;
}

void Solver::setInitialConditions(const Input &input)
{
    using namespace std;
    using namespace boost::property_tree;

    if (input.initialConditionInput().get<std::string>("InitialConditions.type", "") == "restart")
    {
        restartSolution();
        return;
    }

    for (const auto &child: input.initialConditionInput().get_child("InitialConditions"))
    {
        auto scalarFieldIt = scalarFields_.find(child.first);

        if (scalarFieldIt != scalarFields_.end())
        {
            ScalarFiniteVolumeField &field = *scalarFieldIt->second;

            for (const auto &ic: child.second)
            {
                const auto &icTree = ic.second;
                std::string type = icTree.get<string>("type");

                if (type == "circle")
                {
                    Circle circle = Circle(Vector2D(icTree.get<string>("center")), icTree.get<Scalar>("radius"));
                    setCircle(circle, icTree.get<Scalar>("value"), field);
                }
                else if (type == "circleSector")
                {
                    Circle circle = Circle(Vector2D(icTree.get<string>("center")), icTree.get<Scalar>("radius"));
                    Scalar thetaMin = icTree.get<Scalar>("thetaMin");
                    Scalar thetaMax = icTree.get<Scalar>("thetaMax");
                    Scalar innerValue = icTree.get<Scalar>("value");
                    setCircleSector(circle, thetaMin, thetaMax, innerValue, field);
                }
                else if (type == "box")
                {
                    Point2D center = Point2D(icTree.get<string>("center"));
                    Scalar w = icTree.get<Scalar>("width") / 2;
                    Scalar h = icTree.get<Scalar>("height") / 2;

                    Polygon pgn = {
                            Point2D(center.x - w, center.y - h),
                            Point2D(center.x + w, center.y - h),
                            Point2D(center.x + w, center.y + h),
                            Point2D(center.x - w, center.y + h)
                    };

                    setBox(pgn, icTree.get<Scalar>("value"), field);
                }
                else if (type == "uniform")
                    field.fillInterior(icTree.get<Scalar>("value"));
                else if (type == "rotating")
                {
                    setRotating(icTree.get<std::string>("function"),
                                icTree.get<Scalar>("amplitude"),
                                Vector2D(icTree.get<std::string>("center")),
                                field);
                }
                else
                    throw Exception("Input", "setInitialConditions",
                                    "invalid initial condition type \"" + type + "\".");

                printf("Set initial condition \"%s\" of type %s on field \"%s\".\n", ic.first.c_str(), type.c_str(),
                       field.name().c_str());
            }

            continue;
        }

        auto vectorFieldIt = vectorFields_.find(child.first);

        if (vectorFieldIt != vectorFields_.end())
        {
            VectorFiniteVolumeField &field = *vectorFieldIt->second;

            for (const auto &ic: child.second)
            {
                const auto &icTree = ic.second;
                std::string type = icTree.get<string>("type");

                if (type == "circle")
                {
                    Circle circle = Circle(Vector2D(icTree.get<string>("center")), icTree.get<Scalar>("radius"));
                    setCircle(circle, Vector2D(icTree.get<string>("value")), field);
                }
                else if (type == "square")
                {
                    Point2D center = Point2D(icTree.get<string>("center"));
                    Scalar w = icTree.get<Scalar>("width") / 2;
                    Scalar h = icTree.get<Scalar>("height") / 2;

                    Polygon pgn = {
                            Point2D(center.x - w, center.y - h),
                            Point2D(center.x + w, center.y - h),
                            Point2D(center.x + w, center.y + h),
                            Point2D(center.x - w, center.y + h)
                    };

                    setBox(pgn, Vector2D(icTree.get<string>("value")), field);
                }
                else if (type == "uniform")
                    field.fillInterior(Vector2D(icTree.get<string>("value")));
                else if (type == "rotating")
                {
                    setRotating(icTree.get<std::string>("xFunction"),
                                icTree.get<std::string>("yFunction"),
                                Vector2D(icTree.get<std::string>("amplitude")),
                                Vector2D(icTree.get<std::string>("center")),
                                field);
                }
                else
                    throw Exception("Input", "setInitialConditions",
                                    "invalid initial condition type \"" + type + "\".");

                printf("Set initial condition \"%s\" of type %s on field \"%s\".\n", ic.first.c_str(), type.c_str(),
                       field.name().c_str());
            }
        }
    }
}

//- Protected methods

void Solver::setCircle(const Circle &circle, Scalar innerValue, ScalarFiniteVolumeField &field)
{
    const Polygon pgn = circle.polygonize(1000);

    for (const Cell &cell: field.grid().localActiveCells())
    {
        Scalar area = 0.;

        for(const Polygon& xc: intersection(pgn, cell.shape()))
            area += xc.area();

        field(cell) = innerValue * area / cell.shape().area();
    }

    grid_->sendMessages(field);
    interpolateFaces(fv::INVERSE_VOLUME, field);
}

void Solver::setCircle(const Circle &circle, const Vector2D &innerValue, VectorFiniteVolumeField &field)
{
    const Polygon pgn = circle.polygonize(1000);

    for (const Cell &cell: field.grid().localActiveCells())
    {
        Scalar area = 0.;

        for(const Polygon& xc: intersection(pgn, cell.shape()))
            area += xc.area();

        field(cell) = innerValue * area / cell.shape().area();
    }

    grid_->sendMessages(field);
    interpolateFaces(fv::INVERSE_VOLUME, field);
}

void Solver::setCircleSector(const Circle &circle, Scalar thetaMin, Scalar thetaMax, Scalar innerValue,
                             ScalarFiniteVolumeField &field)
{
    thetaMin *= M_PI / 180.;
    thetaMax *= M_PI / 180.;

    while (thetaMin > 2. * M_PI)
        thetaMin -= 2. * M_PI;

    while (thetaMin < 0.)
        thetaMin += 2. * M_PI;

    while (thetaMax > 2. * M_PI)
        thetaMax -= 2. * M_PI;

    while (thetaMax < 0.)
        thetaMax += 2. * M_PI;

    std::vector<Point2D> vtx;
    const int nPts = 1000;
    const Scalar dTheta = (thetaMax - thetaMin) / (nPts - 1);

    for (int i = 0; i < nPts; ++i)
    {
        const Scalar theta = thetaMin + i * dTheta;
        const Vector2D rVec(circle.radius() * std::cos(theta), circle.radius() * std::sin(theta));

        vtx.push_back(circle.centroid() + rVec);
    }

    Polygon pgn(vtx.begin(), vtx.end());

    for (const Cell &cell: field.grid().localActiveCells())
    {
        Scalar area = 0.;

        for(const Polygon& xc: intersection(pgn, cell.shape()))
            area += xc.area();

        field(cell) = innerValue * area / cell.shape().area();
    }

    grid_->sendMessages(field);
    interpolateFaces(fv::INVERSE_VOLUME, field);
}

void Solver::setBox(const Polygon &box, Scalar innerValue, ScalarFiniteVolumeField &field)
{
    for (const Cell &cell: field.grid().localActiveCells())
    {
        Scalar area = 0.;

        for(const Polygon& pgn: intersection(box, cell.shape()))
            area += pgn.area();

        field(cell) = innerValue * area / cell.shape().area();
    }

    grid_->sendMessages(field);
    fv::interpolateFaces(fv::INVERSE_VOLUME, field);
}

void Solver::setBox(const Polygon &box, const Vector2D &innerValue, VectorFiniteVolumeField &field)
{
    for (const Cell &cell: field.grid().localActiveCells())
    {
        Scalar area = 0.;

        for(const Polygon& pgn: intersection(box, cell.shape()))
            area += pgn.area();

        field(cell) = innerValue * area / cell.shape().area();
    }

    grid_->sendMessages(field);
    fv::interpolateFaces(fv::INVERSE_VOLUME, field);
}

void Solver::setRotating(const std::string &function, Scalar amplitude, const Vector2D &center,
                         ScalarFiniteVolumeField &field)
{
    std::function<Scalar(Scalar)> func;

    if (function == "sin")
        func = [](Scalar x) { return sin(x); };
    else if (function == "cos")
        func = [](Scalar x) { return cos(x); };
    else
        throw Exception("Input", "setRotating", "invalid rotation function.");

    for (const Cell &cell: field.grid().cells())
    {
        Vector2D rVec = cell.centroid() - center;
        Scalar theta = atan2(rVec.y, rVec.x);

        field[cell.id()] = amplitude * func(theta);
    }

    for (const Face &face: field.grid().interiorFaces())
    {
        Vector2D rVec = face.centroid() - center;
        Scalar theta = atan2(rVec.y, rVec.x);

        field.faces()[face.id()] = amplitude * func(theta);
    }
}

void Solver::setRotating(const std::string &xFunction, const std::string &yFunction, const Vector2D &amplitude,
                         const Vector2D &center, VectorFiniteVolumeField &field)
{
    std::function<Scalar(Scalar)> xFunc, yFunc;

    if (xFunction == "sin")
        xFunc = [](Scalar x) { return sin(x); };
    else if (xFunction == "cos")
        xFunc = [](Scalar x) { return cos(x); };
    else
        throw Exception("Input", "setRotating", "invalid x rotation function.");

    if (yFunction == "sin")
        yFunc = [](Scalar x) { return sin(x); };
    else if (yFunction == "cos")
        yFunc = [](Scalar x) { return cos(x); };
    else
        throw Exception("Input", "setRotating", "invalid y rotation function.");

    for (const Cell &cell: field.grid().cells())
    {
        Vector2D rVec = cell.centroid() - center;
        Scalar theta = atan2(rVec.y, rVec.x);

        field[cell.id()].x = amplitude.x * xFunc(theta);
        field[cell.id()].y = amplitude.y * yFunc(theta);
    }

    for (const Face &face: field.grid().interiorFaces())
    {
        Vector2D rVec = face.centroid() - center;
        Scalar theta = atan2(rVec.y, rVec.x);

        field.faces()[face.id()].x = amplitude.x * xFunc(theta);
        field.faces()[face.id()].y = amplitude.y * yFunc(theta);
    }
}

void Solver::restartSolution()
{
    using namespace std;
    using namespace boost::filesystem;

    std::regex re("[0-9]+\\.[0-9]+");
    Scalar maxTime = 0.;
    path path;

    for (directory_iterator end, dir("./solution"); dir != end; ++dir)
        if (regex_match(dir->path().filename().string(), re))
        {
            std::smatch match;
            std::regex_search(dir->path().filename().string(), match, re);
            Scalar time = std::stod(match.str());

            if (time > maxTime)
            {
                path = dir->path();
                maxTime = time;
            }
        }

    path /= ("Proc" + std::to_string(grid_->comm().rank())) / "Solution.cgns";

    //- Check to make sure the restart can be done
    if(!exists(path))
        throw Exception("Solver", "restartSolution", "no file \"" + path.string() + "\" needed for restart.");

    int fn;
    cg_open(path.c_str(), CG_MODE_READ, &fn);

    std::vector<Scalar> buffer(grid_->cells().size());
    cgsize_t rmin = 1, rmax = buffer.size();

    for (const auto &field: scalarFields_)
    {
        cg_field_read(fn, 1, 1, 1, field.first.c_str(), CGNS_ENUMV(RealDouble), &rmin, &rmax, field.second->data());
    }

    for (const auto &field: vectorFields_)
    {
        cg_field_read(fn, 1, 1, 1, (field.first + "X").c_str(), CGNS_ENUMV(RealDouble), &rmin, &rmax, buffer.data());

        for (int i = 0; i < buffer.size(); ++i)
            (*field.second)[i].x = buffer[i];

        cg_field_read(fn, 1, 1, 1, (field.first + "Y").c_str(), CGNS_ENUMV(RealDouble), &rmin, &rmax, buffer.data());

        for (int i = 0; i < buffer.size(); ++i)
            (*field.second)[i].y = buffer[i];
    }

    cg_close(fn);
}
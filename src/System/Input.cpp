#include <functional>

#include <boost/property_tree/info_parser.hpp>

#include "Input.h"
#include "Solver.h"
#include "Exception.h"

Input::Input(const std::string &caseDirectory, const std::string &outputPath)
    :
      caseDirectory(caseDirectory),
      outputPath(outputPath)
{

}

void Input::parseInputFile()
{
    using namespace boost::property_tree;

    read_info(caseDirectory + "/case.info", caseInput_);
    read_info(caseDirectory + "/boundaries.info", boundaryInput_);
    read_info(caseDirectory + "/initialConditions.info", initialConditionInput_);
}

void Input::setInitialConditions(const Solver &solver) const
{
    using namespace std;
    using namespace boost::property_tree;

    for(const auto& child: initialConditionInput_.get_child("InitialConditions"))
    {
        auto scalarFieldIt = solver.scalarFields().find(child.first);

        if(scalarFieldIt != solver.scalarFields().end())
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
                else if(type == "uniform")
                    field.fillInterior(icTree.get<Scalar>("value"));
                else if(type == "rotating")
                {
                    setRotating(icTree.get<std::string>("function"),
                                icTree.get<Scalar>("amplitude"),
                                Vector2D(icTree.get<std::string>("center")),
                                field);
                }

                printf("Set initial condition \"%s\" of type %s on field \"%s\".\n", ic.first.c_str(), type.c_str(), field.name.c_str());
            }

            continue;
        }

        auto vectorFieldIt = solver.vectorFields().find(child.first);

        if(vectorFieldIt != solver.vectorFields().end())
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

                printf("Set initial condition \"%s\" of type %s on field \"%s\".\n", ic.first.c_str(), type.c_str(), field.name.c_str());
            }
        }
    }
}

//- Private methods

void Input::setCircle(const Circle &circle, Scalar innerValue, ScalarFiniteVolumeField &field) const
{
    for(const Cell& cell: field.grid.cells())
    {
        if(circle.isInside(cell.centroid()))
            field[cell.id()] = innerValue;
    }

    for(const Face& face: field.grid.interiorFaces())
    {
        if(circle.isInside(face.centroid()))
            field.faces()[face.id()] = innerValue;
    }
}

void Input::setCircle(const Circle &circle, const Vector2D &innerValue, VectorFiniteVolumeField &field) const
{
    for(const Cell& cell: field.grid.cells())
    {
        if(circle.isInside(cell.centroid()))
            field[cell.id()] = innerValue;
    }

    for(const Face& face: field.grid.faces())
    {
        if(circle.isInside(face.centroid()))
            field.faces()[face.id()] = innerValue;
    }
}

void Input::setRotating(const std::string &function, Scalar amplitude, const Vector2D &center, ScalarFiniteVolumeField &field) const
{
    std::function<Scalar(Scalar)> func;

    if(function == "sin")
        func = sin;
    else if(function == "cos")
        func = cos;
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

void Input::setRotating(const std::string &xFunction, const std::string &yFunction, const Vector2D &amplitude, const Vector2D &center, VectorFiniteVolumeField &field) const
{
    std::function<Scalar(Scalar)> xFunc, yFunc;

    if(xFunction == "sin")
        xFunc = sin;
    else if(xFunction == "cos")
        xFunc = cos;
    else
        throw Exception("Input", "setRotating", "invalid x rotation function.");

    if(yFunction == "sin")
        yFunc = sin;
    else if(yFunction == "cos")
        yFunc = cos;
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

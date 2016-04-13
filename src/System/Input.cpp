#include <boost/property_tree/info_parser.hpp>

#include "Input.h"
#include "FiniteVolumeGrid2D.h"

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

void Input::setInitialConditions(const FiniteVolumeGrid2D &grid) const
{
    using namespace std;
    using namespace boost::property_tree;

    for(const auto& child: initialConditionInput_.get_child("InitialConditions"))
    {
        auto scalarFieldIt = grid.scalarFields().find(child.first);

        if(scalarFieldIt != grid.scalarFields().end())
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
                    field.fill(icTree.get<Scalar>("value"));

                printf("Set initial condition \"%s\" of type %s on field \"%s\".\n", ic.first.c_str(), type.c_str(), field.name.c_str());
            }

            continue;
        }

        auto vectorFieldIt = grid.vectorFields().find(child.first);

        if(vectorFieldIt != grid.vectorFields().end())
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
                    field.fill(Vector2D(icTree.get<string>("value")));

                printf("Set initial condition \"%s\" of type %s on field \"%s\".\n", ic.first.c_str(), type.c_str(), field.name.c_str());
            }
        }
    }
}

//- Private methods

void Input::setCircle(const Circle &circle, Scalar innerValue, ScalarFiniteVolumeField &field) const
{
    for(const Cell& cell: field.grid.cells)
    {
        if(circle.isInside(cell.centroid()))
            field[cell.id()] = innerValue;
    }

    for(const Face& face: field.grid.faces)
    {
        if(circle.isInside(face.centroid()))
            field.faces()[face.id()] = innerValue;
    }
}

void Input::setCircle(const Circle &circle, const Vector2D &innerValue, VectorFiniteVolumeField &field) const
{
    for(const Cell& cell: field.grid.cells)
    {
        if(circle.isInside(cell.centroid()))
            field[cell.id()] = innerValue;
    }

    for(const Face& face: field.grid.faces)
    {
        if(circle.isInside(face.centroid()))
            field.faces()[face.id()] = innerValue;
    }
}

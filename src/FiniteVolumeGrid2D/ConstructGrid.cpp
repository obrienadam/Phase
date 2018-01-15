#include <string>

#include "ConstructGrid.h"
#include "StructuredRectilinearGrid.h"
#include "CgnsUnstructuredGrid.h"
#include "Exception.h"

std::shared_ptr<FiniteVolumeGrid2D> constructGrid(const Input &input)
{
    using namespace std;

    //- Grid must be constructed
    string gridType = input.caseInput().get<string>("Grid.type");
    Point2D origin = input.caseInput().get<string>("Grid.origin", "(0,0)");

    if (gridType == "rectilinear")
    {
        Scalar width = input.caseInput().get<Scalar>("Grid.width");
        Scalar height = input.caseInput().get<Scalar>("Grid.height");
        int nCellsX = input.caseInput().get<int>("Grid.nCellsX");
        int nCellsY = input.caseInput().get<int>("Grid.nCellsY");
        Scalar convertToMeters = input.caseInput().get<Scalar>("Grid.convertToMeters", 1.);

        vector<pair<Scalar, Scalar>> xDimRefinements, yDimRefinements;

        Vector2D tmp = input.caseInput().get<string>("Grid.refineX", "(0,0)");
        xDimRefinements.push_back(make_pair(tmp.x, tmp.y));

        tmp = input.caseInput().get<string>("Grid.refineY", "(0,0)");
        yDimRefinements.push_back(make_pair(tmp.x, tmp.y));

        auto grid = std::make_shared<StructuredRectilinearGrid>(width, height,
                                                                nCellsX, nCellsY,
                                                                convertToMeters,
                                                                xDimRefinements,
                                                                yDimRefinements,
                                                                origin);
        return grid;
    }
    else if (gridType == "cgns")
    {
        auto grid = std::make_shared<CgnsUnstructuredGrid>(input);
        return grid;
    }
    else
    {
        throw Exception("", "constructGrid", "invalid grid type \"" + gridType + "\".");
    }
}

std::shared_ptr<FiniteVolumeGrid2D> constructGrid(const Input &input, std::shared_ptr<Communicator> comm)
{
    using namespace std;

    //- Check if a grid needs to be loaded
    if (input.initialConditionInput().get<std::string>("InitialConditions.type", "") == "restart")
    {
        auto grid = std::make_shared<CgnsUnstructuredGrid>();
        grid->loadPartitionedGrid(comm);
        return grid;
    }

    //- Grid must be constructed
    string gridType = input.caseInput().get<string>("Grid.type");
    Point2D origin = input.caseInput().get<string>("Grid.origin", "(0,0)");

    if (gridType == "rectilinear")
    {
        Scalar width = input.caseInput().get<Scalar>("Grid.width");
        Scalar height = input.caseInput().get<Scalar>("Grid.height");
        int nCellsX = input.caseInput().get<int>("Grid.nCellsX");
        int nCellsY = input.caseInput().get<int>("Grid.nCellsY");
        Scalar convertToMeters = input.caseInput().get<Scalar>("Grid.convertToMeters", 1.);

        vector<pair<Scalar, Scalar>> xDimRefinements, yDimRefinements;

        Vector2D tmp = input.caseInput().get<string>("Grid.refineX", "(0,0)");
        xDimRefinements.push_back(make_pair(tmp.x, tmp.y));

        tmp = input.caseInput().get<string>("Grid.refineY", "(0,0)");
        yDimRefinements.push_back(make_pair(tmp.x, tmp.y));

        auto grid = std::make_shared<StructuredRectilinearGrid>(width, height,
                                                                nCellsX, nCellsY,
                                                                convertToMeters,
                                                                xDimRefinements,
                                                                yDimRefinements,
                                                                origin);

        grid->partition(input, comm);
        return grid;
    }
    else if (gridType == "cgns")
    {
        auto grid = std::make_shared<CgnsUnstructuredGrid>(input);
        grid->partition(input, comm);
        return grid;
    }
    else
    {
        throw Exception("", "constructGrid", "invalid grid type \"" + gridType + "\".");
    }
}

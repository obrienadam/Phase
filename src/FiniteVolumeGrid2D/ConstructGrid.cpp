#include <string>

#include "ConstructGrid.h"
#include "StructuredRectilinearGrid.h"
#include "CgnsUnstructuredGrid.h"
#include "Exception.h"

std::shared_ptr<FiniteVolumeGrid2D> constructGrid(const Input& input)
{
    using namespace std;

    string gridType = input.caseInput().get<string>("Grid.type");

    if(gridType == "rectilinear")
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

        return shared_ptr<FiniteVolumeGrid2D>(new StructuredRectilinearGrid(width, height,
                                                                            nCellsX, nCellsY,
                                                                            convertToMeters,
                                                                            xDimRefinements,
                                                                            yDimRefinements));
    }
    else if(gridType == "cgns")
    {
        return shared_ptr<FiniteVolumeGrid2D>(new CgnsUnstructuredGrid(input));
    }
    else
    {
        throw Exception("", "constructGrid", "invalid grid type \"" + gridType + "\".");
    }
}

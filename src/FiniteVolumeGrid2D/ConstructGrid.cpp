#include <string>

#include "ConstructGrid.h"
#include "StructuredRectilinearGrid.h"
#include "CgnsUnstructuredQuadGrid.h"
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

        return shared_ptr<FiniteVolumeGrid2D>(new StructuredRectilinearGrid(width, height, nCellsX, nCellsY, input.caseInput().get<Scalar>("Grid.convertToMeters", 1.)));
    }
    else if(gridType == "cgns")
    {
        return shared_ptr<FiniteVolumeGrid2D>(new CgnsUnstructuredQuadGrid(input));
    }
    else
    {
        throw Exception("", "constructGrid", "invalid grid type \"" + gridType + "\".");
    }
}

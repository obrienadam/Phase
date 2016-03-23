#include <string>

#include "ConstructGrid.h"
#include "RectilinearGrid2D.h"
#include "Exception.h"

std::shared_ptr<FiniteVolumeGrid2D> constructGrid(const Input& input)
{
    using namespace std;

    int nCellsI = input.caseInput().get<int>("Grid.nCellsI");
    int nCellsJ = input.caseInput().get<int>("Grid.nCellsJ");

    string gridType = input.caseInput().get<string>("Grid.type");

    if(gridType == "rectilinear")
    {
        Scalar hx = input.caseInput().get<Scalar>("Grid.spacingX");
        Scalar hy = input.caseInput().get<Scalar>("Grid.spacingY");

        return shared_ptr<RectilinearGrid2D>(new RectilinearGrid2D(nCellsI, nCellsJ, hx, hy));
    }
    else
    {
        throw Exception("", "constructGrid", "invalid grid type \"" + gridType + "\".");
    }
}

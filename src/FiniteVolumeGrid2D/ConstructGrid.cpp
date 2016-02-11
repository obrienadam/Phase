#include <string>

#include "ConstructGrid.h"
#include "EquidistantGrid2D.h"

std::shared_ptr<FiniteVolumeGrid2D> constructGrid(const Input& input)
{
    using namespace std;

    string gridType = input.get<string>("Grid.type");

    if(gridType == "Equidistant")
    {
        int nCellsI = input.get<int>("Grid.nCellsI");
        int nCellsJ = input.get<int>("Grid.nCellsJ");
        Scalar h = input.get<Scalar>("Grid.spacing");

        return shared_ptr<EquidistantGrid2D>(new EquidistantGrid2D(nCellsI, nCellsJ, h));
    }

}

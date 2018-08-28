#include "StructuredGrid2D/StructuredGrid2D.h"

int main()
{
    StructuredGrid2D grid(20, 13, 32, 17);

    for(const Cell& cell: grid.cells())
    {
        std::cout << cell.centroid() << "\n";
    }

    return 0;
}

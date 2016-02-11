#include "DiffusionTerm.h"

DiffusionTerm::DiffusionTerm(const FiniteVolumeGrid2D &grid)
    :
      Term(grid)
{
    Scalar a[5];

    for(int j = 0; j < grid.nCellsJ(); ++j)
        for(int i = 0; i < grid.nCellsI(); ++i)
        {
            int row = grid.cellIndex(i, j);

            if (row == FiniteVolumeGrid2D::INACTIVE)
                continue;

            a[1] = dot(grid.sfe(i, j), grid.sfe(i, j))/dot(grid.rce(i, j), grid.sfe(i, j));
            a[2] = dot(grid.sfw(i, j), grid.sfw(i, j))/dot(grid.rcw(i, j), grid.sfw(i, j));
            a[3] = dot(grid.sfn(i, j), grid.sfn(i, j))/dot(grid.rcn(i, j), grid.sfn(i, j));
            a[4] = dot(grid.sfs(i, j), grid.sfs(i, j))/dot(grid.rcs(i, j), grid.sfs(i, j));
            a[0] = -(a[1] + a[2] + a[3] + a[4]);

            if(grid.inRange(i + 1, j))
                coefficients_.push_back(Triplet(row, grid.cellIndex(i + 1, j), a[1]));

            if(grid.inRange(i - 1, j))
                coefficients_.push_back(Triplet(row, grid.cellIndex(i - 1, j), a[2]));

            if(grid.inRange(i, j + 1))
                coefficients_.push_back(Triplet(row, grid.cellIndex(i, j + 1), a[3]));

            if(grid.inRange(i, j - 1))
                coefficients_.push_back(Triplet(row, grid.cellIndex(i, j - 1), a[4]));

            coefficients_.push_back(Triplet(row, row, a[0]));
        }
}

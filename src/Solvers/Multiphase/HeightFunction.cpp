#include "HeightFunction.h"

HeightFunction::HeightFunction(const Point2D &loc, const Vector2D &normal, Scalar colWidth, Scalar colHeight)
    :
      loc_(loc),
      unitNormal_(normal.unitVec()),
      unitTangent_(unitNormal_.normalVec()),
      colWidth_(colWidth),
      colHeight_(colHeight)
{
    Kernel::Point_2 vertices[4] = {
        (loc_ - colHeight_*unitNormal_/2 - colWidth_*unitTangent_/2).cgalPoint(),
        (loc_ + colHeight_*unitNormal_/2 - colWidth_*unitTangent_/2).cgalPoint(),
        (loc_ + colHeight_*unitNormal_/2 + colWidth_*unitTangent_/2).cgalPoint(),
        (loc_ - colHeight_*unitNormal_/2 + colWidth_*unitTangent_/2).cgalPoint()
    };

    for(int i = 0; i < 4; ++i)
        push_back(vertices[i]);

    init();
}

Scalar HeightFunction::computeHeight(const ScalarFiniteVolumeField &field)
{
    Scalar vof = 0.;

    for(const Cell &cell: field.grid.activeCells())
    {
        const Polygon &pgn = cell.shape();

        if(doIntersect(*this, pgn))
        {
            Scalar overlapArea = intersectionPolygon(*this, pgn).area();
            vof += overlapArea*field[cell.id()];
        }
    }

    return vof/colWidth_;
}

//- External functions

//- External functions
std::vector<HeightFunction> heightFunctionSet(const Point2D& centerLoc, const Vector2D& normal, Scalar colWidth, Scalar colHeight, int nHeightFunctions)
{
    assert(nHeightFunctions > 0 && nHeightFunctions%2 != 0);

    std::vector<HeightFunction> hfs(nHeightFunctions);

    const auto &chf = hfs[nHeightFunctions/2] = HeightFunction(centerLoc, normal, colWidth, colHeight);

    for(int i = 0; i < nHeightFunctions/2; ++i)
    {
        Point2D loc = chf.loc() - (i + 1)*colWidth*chf.tan();
        hfs[nHeightFunctions/2 + 1 + i] = HeightFunction(loc, normal, colWidth, colHeight);

        loc = chf.loc() + (i + 1)*colWidth*chf.tan();
        hfs[nHeightFunctions/2 - 1 - i] = HeightFunction(loc, normal, colWidth, colHeight);
    }

    return hfs;
}

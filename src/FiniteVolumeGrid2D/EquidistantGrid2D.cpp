#include "EquidistantGrid2D.h"
#include "Exception.h"

EquidistantGrid2D::EquidistantGrid2D(int nCellsI, int nCellsJ, Scalar h)
    :
      FiniteVolumeGrid2D(nCellsI, nCellsJ),
      h_(h)
{
    for(int j = 0; j < nNodesJ(); ++j)
        for(int i = 0; i < nNodesI(); ++i)
        {
            cornerNodes_(i, j) = Point2D(i*h_, j*h_);

            if (i < nCellsI && j < nCellsJ)
                cellNodes_(i, j) = Point2D((i + 0.5)*h_, (j + 0.5)*h_);

            if (j < nCellsJ)
                faceNodesI_(i, j) = Point2D(i*h_, (j + 0.5)*h_);

            if (i < nCellsI)
                faceNodesJ_(i, j) = Point2D((i + 0.5)*h_, j*h_);
        }
}

Vector2D EquidistantGrid2D::sf(int i, int j, Face face) const
{   
    switch (face)
    {
    case EAST: return Vector2D(h_, 0.);
    case WEST: return Vector2D(-h_, 0.);
    case NORTH: return Vector2D(0., h_);
    case SOUTH: return Vector2D(0., -h_);
    default: throw Exception("EquidistantGrid2D", "sf", "invalid face.");
    };
}

Vector2D EquidistantGrid2D::rc(int i, int j, Face face) const
{
    switch (face)
    {
    case EAST: return Vector2D(h_, 0.);
    case WEST: return Vector2D(-h_, 0.);
    case NORTH: return Vector2D(0., h_);
    case SOUTH: return Vector2D(0., -h_);
    default: throw Exception("EquidistantGrid2D", "rc", "invalid cell to cell direction.");
    };
}

Vector2D EquidistantGrid2D::rf(int i, int j, Face face) const
{
    switch (face)
    {
    case EAST: return Vector2D(h_/2., 0.);
    case WEST: return Vector2D(-h_/2., 0.);
    case NORTH: return Vector2D(0., h_/2.);
    case SOUTH: return Vector2D(0., -h_/2.);
    default: throw Exception("EquidistantGrid2D", "rf", "invalid cell to face direction.");
    };
}

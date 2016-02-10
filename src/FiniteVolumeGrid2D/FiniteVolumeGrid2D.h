#ifndef FINITE_VOLUME_GRID_2D_H
#define FINITE_VOLUME_GRID_2D_H

#include "Point2D.h"
#include "PointField.h"

class FiniteVolumeGrid2D
{
public:

    enum Face{EAST = 0, WEST = 1, NORTH = 2, SOUTH = 3};
    enum Node{SW = 0, SE = 1, NE = 2, NW = 3};

    FiniteVolumeGrid2D(int nCellsI, int nCellsJ);

    int nCells() const { return cellNodes_.size(); }
    int nCellsI() const { return cellNodes_.sizeI(); }
    int nCellsJ() const { return cellNodes_.sizeJ(); }

    int nNodes() const { return cornerNodes_.size(); }
    int nNodesI() const { return cornerNodes_.sizeI(); }
    int nNodesJ() const { return cornerNodes_.sizeJ(); }

    const Point2D& cellNode(int i, int j) const { return cellNodes_(i, j); }
    const Point2D& cornerNode(int i, int j, Node node) const;
    const Point2D& cornerNode(int i, int j) const { return cornerNodes_(i, j); }

    virtual Vector2D sf(int i, int j, Face face) const = 0;
    virtual Vector2D rf(int i, int j, Face face) const = 0;
    virtual Vector2D rc(int i, int j, Face face) const = 0;

    virtual Scalar cellVolume(int i, int j) const = 0;

    bool inRange(int i, int j) const;

protected:

    PointField cornerNodes_, cellNodes_, faceNodesI_, faceNodesJ_;
};

#endif

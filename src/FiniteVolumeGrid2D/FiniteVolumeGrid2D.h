#ifndef FINITE_VOLUME_GRID_2D_H
#define FINITE_VOLUME_GRID_2D_H

#include "Point2D.h"
#include "PointField.h"
#include "IntegerField.h"

enum Face{EAST = 0, WEST = 1, NORTH = 2, SOUTH = 3};

class FiniteVolumeGrid2D
{
public:

    enum Node{SW = 0, SE = 1, NE = 2, NW = 3};
    enum {INACTIVE = -1};

    FiniteVolumeGrid2D(int nCellsI, int nCellsJ);

    int nCells() const { return cellNodes_.size(); }
    int nCellsI() const { return cellNodes_.sizeI(); }
    int nCellsJ() const { return cellNodes_.sizeJ(); }

    int uCellI() const { return nCellsI() - 1; }
    int uCellJ() const { return nCellsJ() - 1; }

    int nNodes() const { return cornerNodes_.size(); }
    int nNodesI() const { return cornerNodes_.sizeI(); }
    int nNodesJ() const { return cornerNodes_.sizeJ(); }

    const Point2D& cellNode(int i, int j) const { return cellNodes_(i, j); }
    const Point2D& cornerNode(int i, int j, Node node) const;
    const Point2D& cornerNode(int i, int j) const { return cornerNodes_(i, j); }

    virtual Vector2D sf(int i, int j, Face face) const = 0;
    virtual Vector2D rf(int i, int j, Face face) const = 0;
    virtual Vector2D rc(int i, int j, Face face) const = 0;

    Vector2D sfe(int i, int j) const { return sf(i, j, EAST); }
    Vector2D sfw(int i, int j) const { return sf(i, j, WEST); }
    Vector2D sfn(int i, int j) const { return sf(i, j, NORTH); }
    Vector2D sfs(int i, int j) const { return sf(i, j, SOUTH); }

    Vector2D rfe(int i, int j) const { return rf(i, j, EAST); }
    Vector2D rfw(int i, int j) const { return rf(i, j, WEST); }
    Vector2D rfn(int i, int j) const { return rf(i, j, NORTH); }
    Vector2D rfs(int i, int j) const { return rf(i, j, SOUTH); }

    Vector2D rce(int i, int j) const { return rc(i, j, EAST); }
    Vector2D rcw(int i, int j) const { return rc(i, j, WEST); }
    Vector2D rcn(int i, int j) const { return rc(i, j, NORTH); }
    Vector2D rcs(int i, int j) const { return rc(i, j, SOUTH); }

    virtual Scalar cellVolume(int i, int j) const = 0;

    bool inRange(int i, int j) const;
    int cellIndex(int i, int j) const;
    uint nActiveCells() const { return nActiveCells_; }

protected:

    void initCellIndices();

    uint nActiveCells_;

    PointField cornerNodes_, cellNodes_, faceNodesI_, faceNodesJ_;
    IntegerField cellIndices_;
};

#endif

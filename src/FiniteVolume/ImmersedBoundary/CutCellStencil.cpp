#include "CutCellStencil.h"
#include "LineSegment2D.h"

CutCellStencil::CutCellStencil(const Cell &cell, const Shape2D &shape)
        :
        cell_(cell)
{
    auto nodes = cell_.nodes();
    bool intersectedOnce = false;
    Point2D lastXc;

    for (int vtx = 0, endVtx = nodes.size(); vtx < endVtx; ++vtx)
    {
        const Point2D &lNode = nodes[vtx];
        const Point2D &rNode = nodes[(vtx + 1) % nodes.size()];

        bool lNodeIsInside = shape.isInside(lNode);
        bool rNodeIsInside = shape.isInside(rNode);

        LineSegment2D faceGeom(lNode, rNode);

        auto xc = shape.intersections(faceGeom);

        if (!lNodeIsInside && !rNodeIsInside && xc.size() == 0)
            faces_.push_back(LineSegment2D(lNode, rNode));
        else if (!lNodeIsInside && rNodeIsInside && xc.size() == 1)
        {
            faces_.push_back(LineSegment2D(lNode, xc[0]));

            if (intersectedOnce)
                ibFaces_.push_back(LineSegment2D(xc[0], lastXc));
            else
            {
                intersectedOnce = true;
                lastXc = xc[0];
            }
        }
        else if (lNodeIsInside && !rNodeIsInside && xc.size() == 1)
        {
            if (intersectedOnce)
                ibFaces_.push_back(LineSegment2D(lastXc, xc[0]));
            else
            {
                intersectedOnce = true;
                lastXc = xc[0];
            }

            faces_.push_back(LineSegment2D(xc[0], rNode));
        }
        else if (lNodeIsInside && rNodeIsInside && xc.size() == 0)
            continue;
        else
            throw Exception("CutCellStencil", "CutCellStencil", "unexpected truncation case encountered.");
    }
}

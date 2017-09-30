#include "CutCellLink.h"
#include "Exception.h"

CutCellLink::CutCellLink(const Cell &self, const Face &face, const Cell &cell)
        :
        InteriorLink(self, face, cell)
{
    zeta_ = 1.;
    fluidFace_ = LineSegment2D(face.lNode(), face.rNode());
    solidFace_ = LineSegment2D(face.rNode(), face.rNode());
    fluidRFaceVec_ = fluidFace_.center() - self_.centroid();
    solidRFaceVec_ = solidFace_.center() - self_.centroid();
}

CutCellLink::CutCellLink(const Cell &self, const Face &face, const Cell &cell, const Shape2D &shape)
        :
        InteriorLink(self, face, cell)
{
    bool lNodeInside = shape.isInside(face.lNode());
    bool rNodeInside = shape.isInside(face.rNode());

    if (lNodeInside != rNodeInside)
    {
        auto xcs = shape.intersections(LineSegment2D(face.lNode(), face.rNode()));

        if (xcs.size() != 1)
            throw Exception("CutCellLink",
                            "CutCellLink",
                            "Face must intersect the IB exactly once.");

        xc_ = xcs[0];

        if (!lNodeInside)
        {
            fluidFace_ = LineSegment2D(face.lNode(), xc_);
            solidFace_ = LineSegment2D(xc_, face.rNode());
        }
        else
        {
            fluidFace_ = LineSegment2D(xc_, face.rNode());
            solidFace_ = LineSegment2D(face.lNode(), xc_);
        }
    }
    else //- Either both nodes on same side of boundary
    {
        if (!shape.isInside(face.centroid()))
        {
            fluidFace_ = LineSegment2D(face.lNode(), face.rNode());
            solidFace_ = LineSegment2D(face.rNode(), face.rNode());
        }
        else
        {
            fluidFace_ = LineSegment2D(face.lNode(), face.lNode());
            solidFace_ = LineSegment2D(face.lNode(), face.rNode());
        }
    }

    zeta_ = sqrt(fluidFace_.lengthSqr() / outwardNorm_.magSqr());
    fluidRFaceVec_ = fluidFace_.center() - self_.centroid();
    solidRFaceVec_ = solidFace_.center() - self_.centroid();
}

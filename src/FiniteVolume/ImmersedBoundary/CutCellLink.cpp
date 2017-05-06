#include "CutCellLink.h"
#include "Cell.h"
#include "Face.h"
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

    if(lNodeInside != rNodeInside)
    {
        auto xc = shape.intersections(LineSegment2D(face.lNode(), face.rNode()));

        if(xc.size() != 1)
            throw Exception("CutCellLink", "CutCellLink", "face intersects body multiple times.");

        if(!lNodeInside)
        {
            fluidFace_ = LineSegment2D(face.lNode(), xc[0]);
            solidFace_ = LineSegment2D(xc[0], face.rNode());
        }
        else
        {
            fluidFace_ = LineSegment2D(xc[0], face.rNode());
            solidFace_ = LineSegment2D(face.lNode(), xc[0]);
        }
    }
    else //- Either both nodes on same side of boundary
    {
        if(!shape.isInside(face.centroid()))
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

    zeta_ = sqrt(fluidFace_.lengthSqr()/outwardNorm_.magSqr());
    fluidRFaceVec_ = fluidFace_.center() - self_.centroid();
    solidRFaceVec_ = solidFace_.center() - self_.centroid();
}

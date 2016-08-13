#ifndef CELL_H
#define CELL_H

#include <vector>

#include "Types.h"
#include "Polygon.h"

#include "Node.h"
#include "Face.h"
#include "BoundaryFaceLink.h"
#include "InteriorFaceLink.h"
#include "DiagonalCellLink.h"

class Cell
{
public:

    enum {INACTIVE = -1};

    Cell(const std::vector<Label>& nodeIds, const std::vector<Node>& nodes);

    //- Status
    void setActive() { isActive_ = true; }
    void setInactive() { isActive_ = false; globalIndex_ = INACTIVE; }
    bool isActive() const { return isActive_; }

    //- Geometry
    Scalar volume() const { return volume_; }
    const Point2D& centroid() const { return centroid_; }

    //- Id
    Label id() const { return id_; }
    Index globalIndex() const { return globalIndex_; }

    //- Links
    void addDiagonalLink(const Cell& cell);

    const std::vector<InteriorLink>& neighbours() const { return interiorLinks_; }
    const std::vector<BoundaryLink>& boundaries() const { return boundaryLinks_; }
    const std::vector<DiagonalCellLink>& diagonals() const { return diagonalLinks_; }

    //- Nodes
    const std::vector< Ref<const Node>>& nodes() const { return nodes_; }
    const Polygon& shape() const { return cellShape_; }

    Size nFaces() const { return nodes_.size(); }
    Size nInteriorFaces() const { return interiorLinks_.size(); }
    Size nBoundaryFaces() const { return boundaryLinks_.size(); }
    Size nNeighbours() const { return interiorLinks_.size(); }

    bool isInCell(const Point2D& point) const;

private:

    //- Connectivity links, should really only be done by grid classes
    void addBoundaryLink(const Face& face);
    void addInteriorLink(const Face& face, const Cell& cell);

    mutable bool isActive_;
    mutable Index globalIndex_;
    Label id_;

    Polygon cellShape_;

    Scalar volume_;
    Vector2D centroid_;

    std::vector< Ref<const Node> > nodes_;

    std::vector<InteriorLink> interiorLinks_;
    std::vector<BoundaryLink> boundaryLinks_;
    std::vector<DiagonalCellLink> diagonalLinks_;

    friend class FiniteVolumeGrid2D;
};

bool cellsShareFace(const Cell& cellA, const Cell& cellB);

#endif

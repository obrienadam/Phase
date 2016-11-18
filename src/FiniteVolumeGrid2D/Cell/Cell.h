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

class CellGroup;
class FiniteVolumeGrid2D;

class Cell
{
public:

    enum {INACTIVE = -1};

    Cell(const std::vector<Label>& nodeIds, const FiniteVolumeGrid2D& grid);

    //- Status, note the use of mutable types
    void setActive() const { isActive_ = true; }
    void setInactive() const { isActive_ = false; globalIndex_ = INACTIVE; }
    void setFluidCell() const { isFluidCell_ = true; }
    void setNonFluidCell() const { isFluidCell_ = false; }

    bool isActive() const { return isActive_; }
    bool isFluidCell() const { return isFluidCell_; }

    //- Geometry
    Scalar volume() const { return volume_; }
    const Point2D& centroid() const { return centroid_; }

    //- Id
    Label id() const { return id_; }
    void setId(Label id) { id_ = id; }

    Index globalIndex() const { return globalIndex_; }
    void setGlobalIndex(Index index) { globalIndex_ = index; }

    //- Connectivity links, should really only be done by grid classes
    void addDiagonalLink(const Cell& cell);
    void addBoundaryLink(const Face& face);
    void addInteriorLink(const Face& face, const Cell& cell);

    const std::vector<InteriorLink>& neighbours() const { return interiorLinks_; }
    const std::vector<BoundaryLink>& boundaries() const { return boundaryLinks_; }
    const std::vector<DiagonalCellLink>& diagonals() const { return diagonalLinks_; }

    //- Nodes
    const std::vector< Ref<const Node> > nodes() const;
    const Polygon& shape() const { return cellShape_; }

    Size nFaces() const { return cellShape_.vertices().size() - 1; }
    Size nInteriorFaces() const { return interiorLinks_.size(); }
    Size nBoundaryFaces() const { return boundaryLinks_.size(); }
    Size nNeighbours() const { return interiorLinks_.size(); }

    //- Cell group
    const CellGroup& cellGroup() const { return *cellGroupPtr_; }

    bool isInCell(const Point2D& point) const;

private:

    mutable bool isActive_, isFluidCell_; // These flags are just for efficiency
    mutable Index globalIndex_;
    Label id_, globalId_;

    Polygon cellShape_;

    Scalar volume_;
    Vector2D centroid_;

    std::vector<Label> nodeIds_;
    const std::vector<Node>& nodes_;

    std::vector<InteriorLink> interiorLinks_;
    std::vector<BoundaryLink> boundaryLinks_;
    std::vector<DiagonalCellLink> diagonalLinks_;

    const CellGroup* cellGroupPtr_;
};

bool cellsShareFace(const Cell& cellA, const Cell& cellB);

#endif

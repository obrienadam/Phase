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

    enum {INACTIVE = -1, ACTIVE_NO_INDEX = -2};

    Cell(const std::vector<Label>& nodeIds, const FiniteVolumeGrid2D& grid);

    //- Status, note the use of mutable types. Cells can be inactive and still
    //  have a global index defined.
    void setActive() const { localIndex_ = ACTIVE_NO_INDEX; }
    void setInactive() const { localIndex_ = INACTIVE; }
    bool isActive() const { return localIndex_ != INACTIVE; }
    bool isGloballyActive() const { return globalIndex_.size() != 0; }

    //- Geometry
    Scalar volume() const { return volume_; }
    const Point2D& centroid() const { return centroid_; }

    //- Indices
    Index setLocalIndex(Index index) const { return localIndex_ = index; }
    Index localIndex() const { return localIndex_; }

    void setNumberOfGlobalIndices(Size num) const { globalIndex_.resize(num); }
    Index setGlobalIndex(Size num, Index index) const { return globalIndex_[num] = index; }
    Index globalIndex(Size num) const { return globalIndex_[num]; }

    //- Ids
    Label globalId() const { return globalId_; }
    Label setGlobalId(Label id) { return globalId_ = id; }

    Label id() const { return id_; }
    void setId(Label id) { id_ = id; }

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

    mutable Index localIndex_;
    mutable std::vector<Index> globalIndex_; // Indices for linear algebra. May change depending on problem
    Label id_, globalId_; // Indices for identification. Should not normally be changed

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

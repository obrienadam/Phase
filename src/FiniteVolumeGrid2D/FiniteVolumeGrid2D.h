#ifndef FINITE_VOLUME_GRID_2D_H
#define FINITE_VOLUME_GRID_2D_H

#include <vector>

#include "Node.h"
#include "Cell.h"
#include "Face.h"
#include "Patch.h"
#include "BoundingBox.h"

class FiniteVolumeGrid2D
{
public:

    FiniteVolumeGrid2D(size_t nNodes = 0, size_t nCells = 0, size_t nFaces = 0);

    size_t nNodes() const { return nodes_.size(); }
    size_t nCells() const { return cells_.size(); }
    size_t nFaces() const { return faces_.size(); }

    size_t nActiveCells() const { return cells_.size(); }

    size_t createFace(size_t lNodeId, size_t rNodeId, Face::Type type = Face::INTERIOR);
    size_t createCell(const std::vector<size_t>& faceIds);

    bool faceExists(size_t lNodeId, size_t rNodeId) const;
    size_t findFace(size_t lNodeId, size_t rNodeId) const;

    std::string gridInfo() const;

    void computeCellAdjacency();
    size_t computeGlobalIndices();

    void addPatch(const std::string& patchName);
    const std::vector<Patch>& patches() const { return patches_; }

    const Cell& findContainingCell(const Point2D& point, const Cell &guess) const;

    const BoundingBox& boundingBox() const { return bBox_; }

    const std::vector<Cell>& cells() const { return cells_; }
    const std::vector<Node>& nodes() const { return nodes_; }

    const std::vector<Face>& faces() const { return faces_; }
    const std::vector< Ref<const Face> >& interiorFaces() const { return interiorFaces_; }
    const std::vector< Ref<const Face> >& boundaryFaces() const { return boundaryFaces_; }

protected:

    void computeBoundingBox();
    void applyPatch(const std::string& patchName, const std::vector< Ref<Face> >& faces);

    std::vector<Node> nodes_;
    std::vector<Cell> cells_;
    std::vector<Face> faces_;

    std::map<std::pair<size_t, size_t>, size_t> faceDirectory_;

    std::vector< Ref<const Face> > interiorFaces_;
    std::vector< Ref<const Face> > boundaryFaces_;

    std::vector<Patch> patches_;

    BoundingBox bBox_;
};

#endif

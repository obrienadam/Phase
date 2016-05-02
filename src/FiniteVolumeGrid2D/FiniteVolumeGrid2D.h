#ifndef FINITE_VOLUME_GRID_2D_H
#define FINITE_VOLUME_GRID_2D_H

#include <vector>
#include <map>

#include "Node.h"
#include "Cell.h"
#include "Face.h"
#include "Patch.h"
#include "BoundingBox.h"
#include "UniqueCellGroup.h"

class FiniteVolumeGrid2D
{
public:

    FiniteVolumeGrid2D(size_t nNodes = 0, size_t nCells = 0, size_t nFaces = 0);

    //- Size info
    size_t nNodes() const { return nodes_.size(); }
    size_t nCells() const { return cells_.size(); }
    size_t nFaces() const { return faces_.size(); }

    size_t nActiveCells() const { return activeCells_.size(); }

    std::string gridInfo() const;

    //- Create grid entities
    size_t createFace(size_t lNodeId, size_t rNodeId, Face::Type type = Face::INTERIOR);
    size_t createCell(const std::vector<size_t>& faceIds);

    //- Node related methods
    const std::vector<Node>& nodes() const { return nodes_; }

    //- Cell related methods
    const std::vector<Cell>& cells() const { return cells_; }

    const CellGroup& activeCells() const { return activeCells_; }
    const CellGroup& inactiveCells() const { return inactiveCells_; }
    const UniqueCellGroup& fluidCells() const { return fluidCells_; }
    const UniqueCellGroup& cellGroup(const std::string& name) const { return cellGroups_.find(name)->second; }

    const Cell& findContainingCell(const Point2D& point, const Cell &guess) const;

    UniqueCellGroup& moveCellsToFluidCellGroup(const std::vector<size_t>& ids) const;
    UniqueCellGroup& moveAllCellsToFluidCellGroup() const;
    UniqueCellGroup& moveCellsToInactiveCellGroup(const std::vector<size_t>& ids) const;
    UniqueCellGroup& moveCellsToCellGroup(const std::string& name, const std::vector<size_t>& ids) const;

    //- Face related methods
    const std::vector<Face>& faces() const { return faces_; }

    const std::vector< Ref<const Face> >& interiorFaces() const { return interiorFaces_; }
    const std::vector< Ref<const Face> >& boundaryFaces() const { return boundaryFaces_; }

    bool faceExists(size_t lNodeId, size_t rNodeId) const;
    size_t findFace(size_t lNodeId, size_t rNodeId) const { return faceDirectory_.find(std::pair<size_t, size_t>(lNodeId, rNodeId))->second; }

    //- Patch related methods
    void addPatch(const std::string& patchName);
    const std::vector<Patch>& patches() const { return patches_; }

    //- Entity searches
    const Node& findNearestNode(const Point2D& pt) const;

    //- Misc
    const BoundingBox& boundingBox() const { return bBox_; }

protected:

    void initNodes();
    void initCells();

    void constructActiveCellGroup() const;
    void computeBoundingBox();
    void applyPatch(const std::string& patchName, const std::vector< Ref<Face> >& faces);

    //- Node related data
    std::vector<Node> nodes_;

    //- Cell related data
    std::vector<Cell> cells_;

    mutable CellGroup activeCells_; // Stores all currently active cells for which a solution should be computed

    mutable UniqueCellGroup fluidCells_; // Contains all fluid cells for which the standard fv operators will be applied
    mutable UniqueCellGroup inactiveCells_;
    mutable std::map<std::string, UniqueCellGroup> cellGroups_; // Unique cell groups specified by the solver

    //- Face related data
    std::vector<Face> faces_;

    std::map<std::pair<size_t, size_t>, size_t> faceDirectory_; // A directory that can find a face given the two node ids

    std::vector< Ref<const Face> > interiorFaces_; // All interior faces neighboured by two active cells
    std::vector< Ref<const Face> > boundaryFaces_; // All boundary faces neighboured by one active cell

    std::vector<Patch> patches_;

    BoundingBox bBox_;

    //- CGAL stuff
    std::map<Kernel::Point_2, Ref<const Node> > nodeMap_;
};

#endif

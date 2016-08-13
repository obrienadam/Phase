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
#include "Search.h"

class FiniteVolumeGrid2D
{
public:

    FiniteVolumeGrid2D(Size nNodes = 0, Size nCells = 0, Size nFaces = 0);

    //- Size info
    Size nNodes() const { return nodes_.size(); }
    Size nCells() const { return cells_.size(); }
    Size nFaces() const { return faces_.size(); }
    Size nActiveCells() const { return activeCells_.size(); }
    std::string gridInfo() const;

    //- Create grid entities
    Label createCell(const std::vector<Label> &nodeIds);
    Label addNode(const Point2D& point);

    //- Node related methods
    const std::vector<Node>& nodes() const { return nodes_; }

    //- Cell related methods
    const std::vector<Cell>& cells() const { return cells_; }
    const std::vector< Ref<const Cell> >& quadCells() const { return quadCells_; }
    const std::vector< Ref<const Cell> >& triCells() const { return triCells_; }

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

    bool faceExists(Label n1, Label n2) const;
    Label findFace(Label n1, Label n2) const;

    //- Patch related methods
    void addPatch(const std::string& patchName);
    const std::map<std::string, Patch>& patches() const { return patches_; }

    //- Entity searches
    const Node& findNearestNode(const Point2D& pt) const;

    //- Misc
    const BoundingBox& boundingBox() const { return bBox_; }

protected:

    void initNodes();
    void initCells();
    void initConnectivity();

    void constructActiveCellGroup() const;
    void computeBoundingBox();
    void applyPatch(const std::string& patchName, const std::vector< Ref<Face> >& faces);

    //- Node related data
    std::vector<Node> nodes_;

    //- Cell related data
    std::vector<Cell> cells_;

    std::vector< Ref<const Cell> > quadCells_;
    std::vector< Ref<const Cell> > triCells_;

    mutable CellGroup activeCells_; // Stores all currently active cells for which a solution should be computed

    mutable UniqueCellGroup fluidCells_; // Contains all fluid cells for which the standard fv operators will be applied
    mutable UniqueCellGroup inactiveCells_;
    mutable std::map<std::string, UniqueCellGroup> cellGroups_; // Unique cell groups specified by the solver

    //- Face related data
    std::vector<Face> faces_;

    std::map<std::pair<Label, Label>, Label> faceDirectory_; // A directory that can find a face given the two node ids

    std::vector< Ref<const Face> > interiorFaces_; // All interior faces neighboured by two active cells
    std::vector< Ref<const Face> > boundaryFaces_; // All boundary faces neighboured by one active cell

    std::map<std::string, Patch> patches_;

    BoundingBox bBox_;

    //- For node searches
    Search<Node> nodeSearch_;
};

#endif

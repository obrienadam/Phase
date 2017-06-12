#ifndef FINITE_VOLUME_GRID_2D_H
#define FINITE_VOLUME_GRID_2D_H

#include <vector>
#include <map>
#include <unordered_map>

#include "Node.h"
#include "NodeGroup.h"
#include "Cell.h"
#include "CellZone.h"
#include "Face.h"
#include "Patch.h"
#include "BoundingBox.h"
#include "Search.h"
#include "Communicator.h"
#include "Input.h"

class FiniteVolumeGrid2D
{
public:

    FiniteVolumeGrid2D(Size nNodes = 0, Size nCells = 0, Size nFaces = 0);

    //- Initialization
    void init(const std::vector<Point2D> &nodes, const std::vector<Label> &cellInds, const std::vector<Label> &cells);

    void reset();

    //- Size info
    Size nNodes() const
    { return nodes_.size(); }

    Size nCells() const
    { return cells_.size(); }

    Size nFaces() const
    { return faces_.size(); }

    Size nLocalActiveCells() const
    { return localActiveCells().size(); }

    Size nActiveCellsGlobal() const
    { return nActiveCellsGlobal_; }

    std::string gridInfo() const;

    //- Create grid entities
    Label createCell(const std::vector<Label> &nodeIds);

    Label addNode(const Point2D &point);

    //- Node related methods
    const std::vector<Node> &nodes() const
    { return nodes_; }

    const NodeGroup& interiorNodes() const
    { return interiorNodes_; }

    const NodeGroup& boundaryNodes() const
    { return boundaryNodes_; }

    std::vector<Point2D> coords() const;

    std::vector<Scalar> xCoords() const;

    std::vector<Scalar> yCoords() const;

    std::vector<Ref<Node> > getNodes(const std::vector<Label> &ids);

    std::vector<Ref<const Node> > getNodes(const std::vector<Label> &ids) const;

    void assignNodeIds();

    //- Cell related methods

    std::vector<Cell> &cells()
    { return cells_; }

    const std::vector<Cell> &cells() const
    { return cells_; }

    std::vector<int> elementList() const;

    //- Cell groups and zones
    void setCellActive(Cell &cell);

    void setCellsActive(const std::vector<Label> &ids);

    void setCellsActive(const CellGroup &cellGroup);

    void setCellsInactive(const std::vector<Label> &ids);

    void setCellsInactive(const CellGroup &cellGroup);

    void setCellsLocallyInactive(const std::vector<Label> &ids);

    CellGroup &cellGroup(const std::string &name)
    { return cellGroups_.find(name)->second; }

    const CellGroup &cellGroup(const std::string &name) const
    { return cellGroups_.find(name)->second; }

    CellZone &cellZone(const std::string &name)
    { return cellZones_.find(name)->second; }

    const CellZone &cellZone(const std::string &name) const
    { return cellZones_.find(name)->second; }

    CellGroup &createCellGroup(const std::string &name, const std::vector<Label> &ids = std::vector<Label>());

    CellZone &createCellZone(const std::string &name, const std::vector<Label> &ids = std::vector<Label>());

    const std::shared_ptr<CellZone::ZoneRegistry>& cellZoneRegistry() const
    { return cellZoneRegistry_; }

    const CellGroup &localActiveCells() const
    { return localActiveCells_; }

    const CellGroup &globalActiveCells() const
    { return globalActiveCells_; }

    const CellGroup &inactiveCells() const
    { return inactiveCells_; }

    const std::vector<std::shared_ptr<CellZone>> &bufferZones() const
    { return bufferCellZones_; }

    void setAllCellsActive();

    //- Misc cell methods
    const Cell &findContainingCell(const Point2D &point, const Cell &guess) const;

    std::vector<Ref<Cell> > getCells(const std::vector<Label> &ids);

    std::vector<Ref<const Cell> > getCells(const std::vector<Label> &ids) const;

    template<class T>
    std::vector<Label> getCellIds(const T &cells);

    void assignCellIds();

    //- Face related methods

    std::vector<Face> &faces()
    { return faces_; }

    const std::vector<Face> &faces() const
    { return faces_; }

    const FaceGroup& interiorFaces() const
    { return interiorFaces_; }

    const FaceGroup& boundaryFaces() const
    { return boundaryFaces_; }

    bool faceExists(Label n1, Label n2) const;

    Label findFace(Label n1, Label n2) const;

    void assignFaceIds();

    //- Patch related methods
    FaceGroup& createFaceGroup(const std::string& name, const std::vector<Label>& ids = std::vector<Label>());

    Patch& createPatch(const std::string &name, const std::vector<Label> &faces);

    Patch& createPatchByNodes(const std::string &name, const std::vector<Label> &nodes);

    const std::vector<Ref<const Patch>> patches() const;

    FaceGroup& faceGroup(const std::string& name)
    { return faceGroups_.find(name)->second; }

    const FaceGroup& faceGroup(const std::string& name) const
    { return faceGroups_.find(name)->second; }

    Patch& patch(const std::string& name)
    { return patches_.find(name)->second; }

    const Patch& patch(const std::string& name) const
    { return patches_.find(name)->second; }

    const Patch& patch(const Face& face) const
    { return patchRegistry_->find(face.id())->second; }

    //- Entity searches
    const Node &findNearestNode(const Point2D &pt) const;

    std::vector<Ref<Node>> nodesInShape(const Shape2D &shape);

    std::vector<Ref<const Node>> nodesInShape(const Shape2D &shape) const;

    std::vector<std::vector<Ref<const Cell> > > constructSmoothingKernels(Scalar width) const;

    //- Parallel/paritioning
    std::pair<std::vector<int>, std::vector<int>> nodeElementConnectivity() const;

    void partition(const Input &input, const Communicator &comm);

    template<typename T>
    void sendMessages(const Communicator &comm, std::vector<T> &data) const;

    void addNeighbouringProc(int procNo,
                             const std::vector<Label> &sendOrder,
                             const std::vector<Label> &recvOrder);

    //- Active cell ordering, required for lineary algebra!
    void computeGlobalOrdering(const Communicator &comm);

    //- Misc
    const BoundingBox &boundingBox() const
    { return bBox_; }

protected:

    void initNodes();

    void initCells();

    void initConnectivity();

    void computeBoundingBox();

    //- Node related data
    std::vector<Node> nodes_;
    NodeGroup interiorNodes_;
    NodeGroup boundaryNodes_;

    //- Cell related data
    std::vector<Cell> cells_;
    Size nActiveCellsGlobal_;

    //- Default cell groups, active and inactive (import for ordering computations!)
    CellGroup localActiveCells_, inactiveCells_, globalActiveCells_;

    //- Defined cell groups/zones
    std::shared_ptr<CellZone::ZoneRegistry> cellZoneRegistry_;
    std::unordered_map<std::string, CellGroup> cellGroups_;
    std::unordered_map<std::string, CellZone> cellZones_;

    std::vector<int> neighbouringProcs_;
    std::vector<std::shared_ptr<CellGroup>> sendCellGroups_; // shared pointers are used so that zones can be moveable!
    std::vector<std::shared_ptr<CellZone>> bufferCellZones_;

    //- Face related data
    std::vector<Face> faces_;
    std::map<std::pair<Label, Label>, Label> faceDirectory_; // A directory that can find a face given the two node ids

    //- Interior and boundary face data structures
    FaceGroup interiorFaces_;
    FaceGroup boundaryFaces_;

    //- Face groups and patches
    std::shared_ptr<Patch::PatchRegistry> patchRegistry_;
    std::unordered_map<std::string, FaceGroup> faceGroups_;
    std::unordered_map<std::string, Patch> patches_;

    BoundingBox bBox_;
    Polygon domain_;

    //- For node searches
    Search nodeSearch_;
};

#include "FiniteVolumeGrid2D.tpp"

#endif

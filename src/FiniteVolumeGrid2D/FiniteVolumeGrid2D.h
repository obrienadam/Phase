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

    FiniteVolumeGrid2D();

    FiniteVolumeGrid2D(const std::vector<Point2D> &nodes,
                       const std::vector<Label> &cellInds,
                       const std::vector<Label> &cells);

    //- Initialization
    void init(const std::vector<Point2D> &nodes,
              const std::vector<Label> &cellInds,
              const std::vector<Label> &cells);

    void reset();

    //- Size info
    Size nNodes() const
    { return nodes_.size(); }

    Size nCells() const
    { return cells_.size(); }

    Size nFaces() const
    { return faces_.size(); }

    Size nLocalActiveCells() const
    { return localActiveCells_.size(); }

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
    void setCellActive(const Cell &cell);

    template <class const_iterator>
    void setCellsActive(const_iterator begin, const_iterator end)
    {
        localActiveCells_.add(begin, end);
        globalActiveCells_.add(begin, end);
    }

    void setCellInactive(const Cell &cell);

    template <class const_iterator>
    void setCellsInactive(const_iterator begin, const_iterator end)
    {
        localInactiveCells_.add(begin, end);
        globalInactiveCells_.add(begin, end);
    }

    CellGroup &cellGroup(const std::string &name)
    { return cellGroups_.find(name)->second; }

    const CellGroup &cellGroup(const std::string &name) const
    { return cellGroups_.find(name)->second; }

    CellZone &cellZone(const std::string &name)
    { return cellZones_.find(name)->second; }

    const CellZone &cellZone(const std::string &name) const
    { return cellZones_.find(name)->second; }

    CellGroup &createCellGroup(const std::string& name);

    CellZone &createCellZone(const std::string& name, std::shared_ptr<CellZone::ZoneRegistry> registry = nullptr);

    const std::shared_ptr<CellZone::ZoneRegistry>& cellZoneRegistry() const
    { return cellZoneRegistry_; }

    const CellZone &localActiveCells() const
    { return localActiveCells_; }

    const CellZone &globalActiveCells() const
    { return globalActiveCells_; }

    const CellZone &localInactiveCells() const
    { return localActiveCells_; }

    const CellZone &globalInactiveCells() const
    { return globalInactiveCells_; }

    const std::vector<CellZone> &bufferZones() const
    { return bufferCellZones_; }

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

    bool hasPatch(const std::string& name) const
    { return patches_.find(name) != patches_.end(); }

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

    //- Parallel/paritioning
    const Communicator &comm() const
    { return *comm_; }

    std::pair<std::vector<int>, std::vector<int>> nodeElementConnectivity() const;

    void partition(const Input &input, std::shared_ptr<Communicator> comm);

    template<class T>
    void sendMessages(std::vector<T> &data) const;

    //- Active cell ordering, required for lineary algebra!
    void computeGlobalOrdering();

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

    //- Local cell zones
    CellZone localActiveCells_, localInactiveCells_;

    //- Cells on this proc that are globally active/inactive
    CellZone globalActiveCells_, globalInactiveCells_;

    //- User defined cell groups/zones
    std::shared_ptr<CellZone::ZoneRegistry> cellZoneRegistry_;
    std::unordered_map<std::string, CellGroup> cellGroups_;
    std::unordered_map<std::string, CellZone> cellZones_;

    //- Communication zones
    std::shared_ptr<Communicator> comm_;
    std::vector<CellGroup> sendCellGroups_; // shared pointers are used so that zones can be moveable!
    std::vector<CellZone> bufferCellZones_;

    //- Face related data
    std::vector<Face> faces_;
    std::map<std::pair<Label, Label>, Label> faceDirectory_; // A directory that can find a face given the two node ids

    //- Interior and boundary face data structures
    FaceGroup interiorFaces_;
    FaceGroup boundaryFaces_;

    //- User defined face groups and patches
    std::shared_ptr<Patch::PatchRegistry> patchRegistry_;
    std::unordered_map<std::string, FaceGroup> faceGroups_;
    std::unordered_map<std::string, Patch> patches_;

    BoundingBox bBox_;

    //- For node searches
    Search nodeSearch_;
};

#include "FiniteVolumeGrid2D.tpp"

#endif

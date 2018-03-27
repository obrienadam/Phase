#ifndef FINITE_VOLUME_GRID_2D_H
#define FINITE_VOLUME_GRID_2D_H

#include "System/Input.h"
#include "System/Communicator.h"

#include "Node/Node.h"
#include "Node/NodeGroup.h"
#include "Cell/Cell.h"
#include "Cell/CellZone.h"
#include "Face/Face.h"
#include "Face/Patch.h"

#include "Geometry/BoundingBox.h"

class FiniteVolumeGrid2D
{
public:

    FiniteVolumeGrid2D();

    FiniteVolumeGrid2D(const std::vector<Point2D> &nodes,
                       const std::vector<Label> &cellInds,
                       const std::vector<Label> &cells,
                       const Point2D &origin);

    //- Initialization
    virtual void init(const std::vector<Point2D> &nodes,
                      const std::vector<Label> &cellInds,
                      const std::vector<Label> &cells,
                      const Point2D &origin);

    virtual void reset();

    //- Size info
    Size nNodes() const
    { return nodes_.size(); }

    Size nCells() const
    { return cells_.size(); }

    std::string gridInfo() const;

    //- Create grid entities
    virtual Label createCell(const std::vector<Label> &nodeIds);

    virtual Label addNode(const Point2D &point);

    //- Node related methods
    const std::vector<Node> &nodes() const
    { return nodes_; }

    const NodeGroup &interiorNodes() const
    { return interiorNodes_; }

    const NodeGroup &boundaryNodes() const
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
    template<class const_iterator>
    void setCellsActive(const_iterator begin, const_iterator end)
    {
        localActiveCells_.add(begin, end);
    }

    template<class const_iterator>
    void setCellsInActive(const_iterator begin, const_iterator end)
    {
        localInactiveCells_.add(begin, end);
    }

    void setCellActive(const Cell &cell);

    void setCellInactive(const Cell &cell);

    void updateGlobalActiveCells();

    template<class const_iterator>
    void setCellsInactive(const_iterator begin, const_iterator end)
    {
        localInactiveCells_.add(begin, end);
        globalInactiveCells_.add(begin, end);
    }

    //- Cell group access
    CellGroup &cellGroup(const std::string &name)
    { return *cellGroups_.find(name)->second; }

    const CellGroup &cellGroup(const std::string &name) const
    { return *cellGroups_.find(name)->second; }

    std::shared_ptr<CellGroup> &cellGroupPtr(const std::string &name)
    { return cellGroups_.find(name)->second; }

    //- Cell zone access
    CellZone &cellZone(const std::string &name)
    { return *cellZones_.find(name)->second; }

    const CellZone &cellZone(const std::string &name) const
    { return *cellZones_.find(name)->second; }

    std::shared_ptr<CellZone> &cellZonePtr(const std::string &name)
    { return cellZones_.find(name)->second; }

    CellGroup globalCellGroup(const CellGroup &localGroup) const;

    //- Cell group creation
    CellGroup &createCellGroup(const std::string &name);

    //- Cell zone creation
    CellZone &createCellZone(const std::string &name, std::shared_ptr<CellZone::ZoneRegistry> registry = nullptr);

    const std::shared_ptr<CellZone::ZoneRegistry> &cellZoneRegistry() const
    { return cellZoneRegistry_; }

    const CellZone &localActiveCells() const
    { return localActiveCells_; }

    const CellZone &localInactiveCells() const
    { return localActiveCells_; }

    const CellZone &globalActiveCells() const
    { return globalActiveCells_; }

    const CellZone &globalInactiveCells() const
    { return globalInactiveCells_; }

    const std::vector<CellZone> &bufferZones() const
    { return bufferCellZones_; }

    //- Face related methods

    std::vector<Face> &faces()
    { return faces_; }

    const std::vector<Face> &faces() const
    { return faces_; }

    const FaceGroup &interiorFaces() const
    { return interiorFaces_; }

    const FaceGroup &boundaryFaces() const
    { return boundaryFaces_; }

    bool faceExists(Label n1, Label n2) const;

    Label findFace(Label n1, Label n2) const;

    void assignFaceIds();

    //- Patch related methods
    FaceGroup &createFaceGroup(const std::string &name, const std::vector<Label> &ids = std::vector<Label>());

    Patch &createPatch(const std::string &name, const std::vector<Label> &faces);

    Patch &createPatchByNodes(const std::string &name, const std::vector<Label> &nodes);

    const std::vector<Ref<const Patch>> patches() const;

    FaceGroup &faceGroup(const std::string &name)
    { return faceGroups_.find(name)->second; }

    const FaceGroup &faceGroup(const std::string &name) const
    { return faceGroups_.find(name)->second; }

    bool hasPatch(const std::string &name) const
    { return patches_.find(name) != patches_.end(); }

    Patch &patch(const std::string &name)
    { return patches_.find(name)->second; }

    const Patch &patch(const std::string &name) const
    { return patches_.find(name)->second; }

    const Patch &patch(const Face &face) const
    { return patchRegistry_->find(face.id())->second; }

    //- Entity searches
    const Node &findNearestNode(const Point2D &pt) const;

    std::vector<Ref<const Node>> findNearestNodes(const Point2D &pt, int nNodes) const;

    //- Parallel/paritioning
    const Communicator &comm() const
    { return *comm_; }

    std::pair<std::vector<int>, std::vector<int>> nodeElementConnectivity() const;

    std::unordered_map<std::string, std::vector<int>> patchToNodeMap() const;

    void partition(const Input &input, std::shared_ptr<Communicator> comm);

    template<class T>
    void sendMessages(std::vector<T> &data) const;

    template<class T>
    void sendMessages(std::vector<T> &data, Size nSets) const;

    //- Misc
    const BoundingBox &boundingBox() const
    { return bBox_; }

protected:

    void initConnectivity();

    void computeBoundingBox();

    //- Node related data
    std::vector<Node> nodes_;
    NodeGroup interiorNodes_, boundaryNodes_, nodeGroup_;

    //- Cell related data
    std::vector<Cell> cells_;

    //- Local cell zones
    CellZone localActiveCells_, localInactiveCells_;

    //- Cells on this proc that are globally active/inactive
    CellZone globalActiveCells_, globalInactiveCells_;

    //- User defined cell groups/zones
    std::shared_ptr<CellZone::ZoneRegistry> cellZoneRegistry_;
    std::unordered_map<std::string, std::shared_ptr<CellGroup>> cellGroups_;
    std::unordered_map<std::string, std::shared_ptr<CellZone>> cellZones_;

    //- Communication zones
    std::shared_ptr<Communicator> comm_;
    std::vector<CellGroup> sendCellGroups_;
    std::vector<CellZone> bufferCellZones_;

    //- Face related data
    std::vector<Face> faces_;
    std::map<std::pair<Label, Label>, Label> faceDirectory_; // A directory that can find a face given the two node ids

    //- Interior and boundary face data structures
    FaceGroup interiorFaces_, boundaryFaces_;

    //- User defined face groups and patches
    std::shared_ptr<Patch::PatchRegistry> patchRegistry_;
    std::unordered_map<std::string, FaceGroup> faceGroups_;
    std::unordered_map<std::string, Patch> patches_;

    BoundingBox bBox_;
};

#include "FiniteVolumeGrid2D.tpp"

#endif

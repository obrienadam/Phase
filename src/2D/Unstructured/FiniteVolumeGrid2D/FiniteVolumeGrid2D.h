#ifndef PHASE_FINITE_VOLUME_GRID_2D_H
#define PHASE_FINITE_VOLUME_GRID_2D_H

#include <unordered_map>

#include "System/Communicator.h"
#include "System/Input.h"

#include "Cell/Cell.h"
#include "Cell/CellGroup.h"
#include "Face/Face.h"
#include "Face/FaceGroup.h"
#include "Node/Node.h"
#include "Node/NodeGroup.h"

#include "Geometry/BoundingBox.h"

class FiniteVolumeGrid2D {
public:
  FiniteVolumeGrid2D();

  FiniteVolumeGrid2D(const std::vector<Point2D> &nodes,
                     const std::vector<Label> &cptr,
                     const std::vector<Label> &cind, const Point2D &origin);

  //- Make non-copyable
  FiniteVolumeGrid2D(const FiniteVolumeGrid2D &grid) = delete;

  FiniteVolumeGrid2D &operator=(const FiniteVolumeGrid2D &grid) = delete;

  //- Initialization
  virtual void init(const std::vector<Point2D> &nodes,
                    const std::vector<Label> &cptr,
                    const std::vector<Label> &cind, const Point2D &origin);

  virtual void reset();

  //- Size info
  Size nNodes() const { return nodes_.size(); }

  Size nCells() const { return cells_.size(); }

  Size nFaces() const { return faces_.size(); }

  std::string info() const;

  //- Create grid entities
  Label createCell(const std::vector<Label> &nodeIds);

  Label addNode(const Point2D &point);

  //- Node related methods
  const std::vector<Node> &nodes() const { return nodes_; }

  std::vector<Point2D> coords() const {
    return std::vector<Point2D>(nodes_.begin(), nodes_.end());
  }

  //- Cell related methods
  const std::vector<Cell> &cells() const { return cells_; }

  const CellGroup &localCells() const { return localCells_; }

  const CellGroup &globalCells() const { return globalCells_; }

  const std::vector<Label> &cellOwnership() const { return cellOwnership_; }

  const std::vector<Label> &globalIds() const { return globalIds_; }

  CellGroup globalCellGroup(const CellGroup &localGroup) const;

  const std::vector<CellGroup> &sendGroups() const { return sendCellGroups_; }

  const std::vector<CellGroup> &bufferGroups() const {
    return bufferCellGroups_;
  }

  //- Face related methods
  std::vector<Face> &faces() { return faces_; }

  const std::vector<Face> &faces() const { return faces_; }

  const FaceGroup &interiorFaces() const { return interiorFaces_; }

  const FaceGroup &boundaryFaces() const { return boundaryFaces_; }

  bool faceExists(Label n1, Label n2) const;

  Label findFace(Label n1, Label n2) const;

  //- Patch related methods
  FaceGroup &createPatch(const std::string &name,
                         const std::vector<Label> &faces);

  FaceGroup &createPatchByNodes(const std::string &name,
                                const std::vector<Label> &nodes);

  std::vector<Ref<const FaceGroup>> patches() const;

  FaceGroup &patch(const std::string &name) {
    return patches_.find(name)->second;
  }

  const FaceGroup &patch(const std::string &name) const {
    return patches_.find(name)->second;
  }

  const FaceGroup &patch(const Face &face) const {
    return patchRegistry_.find(face.id())->second;
  }

  //- Entity searches
  const Node &findNearestNode(const Point2D &pt) const;

  std::vector<Ref<const Node>> findNearestNodes(const Point2D &pt,
                                                int nNodes) const;

  //- Parallel/paritioning
  const Communicator &comm() const { return *comm_; }

  std::vector<int> eptr() const;

  std::vector<int> eind() const;

  std::pair<std::vector<int>, std::vector<int>> connectivityGraph() const;

  std::unordered_map<std::string, std::vector<int>> patchToNodeMap() const;

  void partition(const Input &input);

  template <class T> void sendMessages(std::vector<T> &data) const;

  template <class T> void sendMessages(std::vector<T> &data, Size nSets) const;

  //- Misc
  const BoundingBox &boundingBox() const { return bBox_; }

protected:
  void init();

  void initPatches(
      const std::unordered_map<std::string, std::vector<Label>> &patches);

  void initCommBuffers(const std::vector<Label> &ownership,
                       const std::vector<Label> &globalIds);

  //- Node related data
  std::vector<Node> nodes_;

  NodeGroup nodeGroup_;

  //- Cell related data
  std::vector<Cell> cells_;

  //- Cell groups for searching
  CellGroup localCells_;

  CellGroup globalCells_;

  std::vector<Label> cellOwnership_, globalIds_;

  //- Communication zones
  std::shared_ptr<Communicator> comm_;

  std::vector<CellGroup> sendCellGroups_, bufferCellGroups_;

  //- Face related data
  std::vector<Face> faces_;

  std::map<std::pair<Label, Label>, Label>
      faceDirectory_; // A directory that can find a face given the two node ids

  //- Interior and boundary face data structures
  FaceGroup interiorFaces_, boundaryFaces_;

  std::unordered_map<std::string, FaceGroup> patches_;

  std::unordered_map<Label, Ref<const FaceGroup>> patchRegistry_;

  BoundingBox bBox_;
};

#include "FiniteVolumeGrid2D.tpp"

#endif

#include <numeric>

#include <metis.h>

#include "FiniteVolumeGrid2D.h"

FiniteVolumeGrid2D::FiniteVolumeGrid2D()
    :
      interiorFaces_("InteriorFaces"),
      boundaryFaces_("BoundaryFaces"),
      localCells_("LocalCells"),
      globalCells_("GlobalCells"),
      comm_(std::make_shared<Communicator>())
{

}

FiniteVolumeGrid2D::FiniteVolumeGrid2D(const std::vector<Point2D> &nodes,
                                       const std::vector<Label> &cptr,
                                       const std::vector<Label> &cind,
                                       const Point2D &origin)
    :
      FiniteVolumeGrid2D()
{
    init(nodes, cptr, cind, origin);
}

void FiniteVolumeGrid2D::init(const std::vector<Point2D> &nodes,
                              const std::vector<Label> &cptr,
                              const std::vector<Label> &cind,
                              const Point2D &origin)
{
    reset();

    for (const Point2D &node: nodes)
        nodes_.push_back(Node(node + origin, *this));

    cells_.reserve(cptr.size() - 1); // very important, can break without reserve
    for (int i = 0; i < cptr.size() - 1; ++i)
        createCell(std::vector<Label>(cind.begin() + cptr[i], cind.begin() + cptr[i + 1]));

    init();
}

void FiniteVolumeGrid2D::reset()
{
    //- Node related data
    nodes_.clear();
    nodeGroup_.clear();

    //- Cell related data
    cells_.clear();
    localCells_.clear();
    globalCells_.clear();

    //- Communication zones
    sendCellGroups_.clear(); // shared pointers are used so that zones can be moveable!
    bufferCellGroups_.clear();

    //- Face related data
    faces_.clear();
    faceDirectory_.clear(); // A directory that can find a face given the two node ids

    //- Interior and boundary face data structures
    interiorFaces_.clear();
    boundaryFaces_.clear();

    //- User defined face groups and patches
    patches_.clear();
    bBox_ = BoundingBox(Point2D(0., 0.), Point2D(0., 0.));
}

//- size info
std::string FiniteVolumeGrid2D::info() const
{
    using namespace std;

    return "Finite volume grid info:\n"
           "------------------------\n"
           "Number of nodes: " + to_string(nodes_.size()) + "\n" +
            "Number of cells total: " + to_string(cells_.size()) + "\n" +
            "Number of interior faces: " + to_string(interiorFaces_.size()) + "\n" +
            "Number of boundary faces: " + to_string(boundaryFaces_.size()) + "\n" +
            "Number of faces total: " + to_string(faces_.size()) + "\n";
}

//- Create grid entities
Label FiniteVolumeGrid2D::createCell(const std::vector<Label> &nodeIds)
{
    using namespace std;

    cells_.push_back(Cell(nodeIds, *this));
    Cell &newCell = cells_.back();

    for (Label id: nodeIds)
        nodes_[id].addCell(newCell);

    for (Label i = 0, end = nodeIds.size(); i < end; ++i)
    {
        Label n1 = nodeIds[i], n2 = nodeIds[(i + 1) % end];

        auto key = n1 < n2 ? make_pair(n1, n2) : make_pair(n2, n1);
        auto it = faceDirectory_.find(key);

        if (it == faceDirectory_.end()) // face doesn't exist, so create it
        {
            faces_.push_back(Face(n1, n2, *this, Face::BOUNDARY));
            Face &face = faces_.back();
            faceDirectory_[key] = face.id();
            face.addCell(newCell);
        }
        else // face already exists, but is now an interior face
        {
            Face &face = faces_[it->second];
            face.setType(Face::INTERIOR);
            face.addCell(newCell);
        }
    }

    return newCell.id();
}

Label FiniteVolumeGrid2D::addNode(const Point2D &point)
{
    nodes_.push_back(Node(point, *this));
    return nodes_.back().id();
}

CellGroup FiniteVolumeGrid2D::globalCellGroup(const CellGroup &localGroup) const
{
    if (!comm_ || comm_->nProcs() == 1)
        return localGroup;

    CellGroup globalGroup(localGroup);
    std::vector<int> isInGlobalGroup(cells_.size());

    std::transform(cells_.begin(), cells_.end(), isInGlobalGroup.begin(), [&globalGroup](const Cell &cell) -> int
    {
        return (int) globalGroup.isInSet(cell);
    });

    sendMessages(isInGlobalGroup);

    for (const CellGroup &buffer: bufferGroups())
        for (const Cell &cell: buffer)
            if (isInGlobalGroup[cell.id()])
                globalGroup.add(cell);

    return globalGroup;
}

//- Face related methods

bool FiniteVolumeGrid2D::faceExists(Label n1, Label n2) const
{
    using namespace std;

    return faceDirectory_.find(n1 < n2 ? make_pair(n1, n2) : make_pair(n2, n1)) != faceDirectory_.end();
}

Label FiniteVolumeGrid2D::findFace(Label n1, Label n2) const
{
    using namespace std;

    auto it = faceDirectory_.find(n1 < n2 ? make_pair(n1, n2) : make_pair(n2, n1));

    if (it == faceDirectory_.end())
        throw Exception("FiniteVolumeGrid2D", "findFace",
                        "no face found between n1 = " + to_string(n1) + ", n2 = " + to_string(n2) + ".");

    return it->second;
}

//- Patch related methods
FaceGroup &FiniteVolumeGrid2D::createPatch(const std::string &name, const std::vector<Label> &faces)
{
    FaceGroup &patch = patches_[name] = FaceGroup(name);

    for(Label fid: faces)
    {
        auto insert = patchRegistry_.insert(std::make_pair(fid, std::cref(patch)));

        if(!insert.second)
            insert.first->second = std::cref(patch);

        patch.add(faces_[fid]);
    }

    return patch;
}

FaceGroup &FiniteVolumeGrid2D::createPatchByNodes(const std::string &name, const std::vector<Label> &nodes)
{
    FaceGroup &patch = patches_[name] = FaceGroup(name);

    for (Label lid = 0, rid = 1; rid < nodes.size(); lid += 2, rid += 2)
    {
        Label fid = findFace(nodes[lid], nodes[rid]);
        patch.add(faces_[fid]);

        auto insert = patchRegistry_.insert(std::make_pair(fid, std::cref(patch)));

        if(!insert.second)
            insert.first->second = std::cref(patch);
    }

    return patch;
}

std::vector<Ref<const FaceGroup> > FiniteVolumeGrid2D::patches() const
{
    std::vector<Ref<const FaceGroup>> patches;
    patches.reserve(patches.size());

    for(const auto &entry: patches_)
        patches.push_back(std::cref(entry.second));

    return patches;
}

const Node &FiniteVolumeGrid2D::findNearestNode(const Point2D &pt) const
{
    return nodeGroup_.nearestItems(pt, 1)[0];
}

std::vector<Ref<const Node>> FiniteVolumeGrid2D::findNearestNodes(const Point2D &pt, int nNodes) const
{
    return nodeGroup_.nearestItems(pt, nNodes);
}

std::vector<int> FiniteVolumeGrid2D::eptr() const
{
    std::vector<int> eptr(1, 0);
    eptr.reserve(cells_.size() + 1);

    std::transform(cells_.begin(), cells_.end(), std::back_inserter(eptr), [](const Cell &cell)
    { return cell.nodes().size(); });

    std::partial_sum(eptr.begin(), eptr.end(), eptr.begin());

    return eptr;
}

std::vector<int> FiniteVolumeGrid2D::eind() const
{
    std::vector<int> eind;
    eind.reserve(4 * cells_.size());

    for (const Cell &cell: cells_)
        std::transform(cell.nodes().begin(), cell.nodes().end(), std::back_inserter(eind),
                       [](const Node &node) { return node.id(); });

    return eind;
}

std::pair<std::vector<int>, std::vector<int>> FiniteVolumeGrid2D::connectivityGraph() const
{
    auto graph = std::make_pair(std::vector<int>(1, 0), std::vector<int>());

    for (const Cell &cell: cells_)
    {
        graph.first.push_back(graph.first.back() + cell.neighbours().size());
        for (const auto &nb: cell.neighbours())
            graph.second.push_back(nb.cell().id());
    }

    return graph;
}

std::unordered_map<std::string, std::vector<int>> FiniteVolumeGrid2D::patchToNodeMap() const
{
    using namespace std;
    unordered_map<string, vector<int>> patchToNodes;

    for (const FaceGroup &patch: patches())
    {
        vector<int> nodes;
        nodes.reserve(2 * patch.items().size());

        for (const Face &face: patch)
        {
            nodes.push_back(face.lNode().id());
            nodes.push_back(face.rNode().id());
        }

        patchToNodes[patch.name()] = nodes;
    }

    return patchToNodes;
}

void FiniteVolumeGrid2D::partition(const Input &input)
{
    using namespace std;

    if (comm_->nProcs() == 1) // no need to perform a partition
        return;

    comm_->printf("Partitioning grid into %d partitions...\n", comm_->nProcs());
    vector<idx_t> cellPartition(nCells());

    if (comm_->isMainProc()) // partition is performed on main proc
    {
        idx_t nPartitions = comm_->nProcs();
        idx_t nElems = nCells();
        idx_t nNodes = this->nNodes();
        idx_t nCommon = 2; //- face connectivity weighting only
        idx_t objVal;
        vector<idx_t> nodePartition(this->nNodes());

        int status = METIS_PartMeshDual(&nElems, &nNodes,
                                        eptr().data(),
                                        eind().data(),
                                        NULL, NULL,
                                        &nCommon, &nPartitions,
                                        NULL, NULL, &objVal,
                                        cellPartition.data(), nodePartition.data());
        if (status == METIS_OK)
            comm_->printf("Sucessfully computed partitioning.\n");
        else
            throw Exception("FiniteVolumeGrid2D", "partition", "an error occurred during partitioning.");
    }

    //- Broadcast the partitioning to other processes
    comm_->broadcast(comm_->mainProcNo(), cellPartition);

    //- Criteria to see if a cell is retained on a particular proc
    auto addCellToThisProc = [this, &cellPartition](const Cell &cell, Scalar r = 0.) -> bool
    {
        if (cellPartition[cell.id()] == comm_->rank())
            return true;

        for (const CellLink &nb: cell.cellLinks())
            if (cellPartition[nb.cell().id()] == comm_->rank())
                return true;

        for (const Cell &kCell: globalCells_.itemsWithin(Circle(cell.centroid(), r)))
            if (cellPartition[kCell.id()] == comm_->rank())
                return true;

        return false;
    };

    //- Construct the crs representation of the local grid
    comm_->printf("Computing the local cell domains...\n");
    vector<Point2D> nodes;
    vector<Label> cellInds(1, 0), cellNodeIds, cellProc;
    unordered_map<Label, Label> cellGlobalToLocalIdMap, cellLocalToGlobalIdMap;
    vector<int> localNodeId(nodes_.size(), -1);
    Scalar r = input.caseInput().get<Scalar>("Grid.minBufferWidth", 0.);

    for (const Cell &cell: cells_)
        if (addCellToThisProc(cell, r))
        {
            cellInds.push_back(cellInds.back() + cell.nodes().size());
            cellProc.push_back(cellPartition[cell.id()]);

            cellGlobalToLocalIdMap[cell.id()] = cellInds.size() - 2;
            cellLocalToGlobalIdMap[cellInds.size() - 2] = cell.id();

            for (const Node &node: cell.nodes())
            {
                if (localNodeId[node.id()] == -1)
                {
                    localNodeId[node.id()] = nodes.size();
                    nodes.push_back(node);
                }

                cellNodeIds.push_back(localNodeId[node.id()]);
            }
        }

    //- Boundary patches
    comm_->printf("Computing the local boundary patches...\n");
    std::unordered_map<std::string, std::vector<Label>> localPatches;

    for (const FaceGroup &patch: patches())
    {
        vector<Label> nodeIds;

        for (const Face &face: patch)
        {
            int lid = localNodeId[face.lNode().id()];
            int rid = localNodeId[face.rNode().id()];

            if (lid > -1 && rid > -1)
                nodeIds.insert(nodeIds.end(), {(Label) lid, (Label) rid});
        }

        if (!nodeIds.empty())
            localPatches[patch.name()] = nodeIds;
    }

    //- Now re-initialize local domains
    comm_->printf("Initializing local domains...\n");

    init(nodes, cellInds, cellNodeIds, Point2D(0., 0.));
    initPatches(localPatches);

    comm_->printf("Finished initializing local domains.\n");

    auto nLocalCells = comm_->allGather(localCells_.size());
    auto gid = std::accumulate(nLocalCells.begin(), nLocalCells.begin() + comm_->rank(), 0);

    for (const Cell &cell: cells_)
    {
        cellOwnership_[cell.id()] = cellPartition[cellLocalToGlobalIdMap[cell.id()]];
        globalIds_[cell.id()] = cellLocalToGlobalIdMap[cell.id()];
    }

    comm_->printf("Initiating inter-process communication buffers...\n");

    initCommBuffers(cellOwnership_, globalIds_);
}

//- Protected methods

void FiniteVolumeGrid2D::init()
{
    nodeGroup_ = NodeGroup(nodes_.begin(), nodes_.end(), "NodeGroup");

    nodeGroup_.clear();
    nodeGroup_.add(nodes_.begin(), nodes_.end());

    boundaryFaces_.clear();
    interiorFaces_.clear();

    for (const Face &face: faces_)
        if (face.isBoundary())
        {
            Cell &cell = cells_[face.lCell().id()];
            cell.addBoundaryLink(face);

            boundaryFaces_.add(face);
        }
        else
        {
            Cell &lCell = cells_[face.lCell().id()];
            Cell &rCell = cells_[face.rCell().id()];

            lCell.addInteriorLink(face, rCell);
            rCell.addInteriorLink(face, lCell);

            interiorFaces_.add(face);
        }

    //- Initialize diagonal links
    for (Cell &cell: cells_)
        for (const Node &node: cell.nodes())
            for (const Cell &kCell: node.cells())
            {
                if (&cell == &kCell)
                    continue;
                else if (!cellsShareFace(cell, kCell))
                    cell.addDiagonalLink(kCell);
            }

    //    //- Initialize the patch registry
    //    patchRegistry_.clear();
    //
    //    for (const auto &entry: patches_)
    //        for (const Face &face: entry.second)
    //            patchRegistry_[face.id()] = std::cref(entry.second);

    //- Init the local and global cell groups
    localCells_.clear();
    localCells_.add(cells_.begin(), cells_.end());

    globalCells_.clear();
    globalCells_.add(localCells_);

    cellOwnership_.resize(globalCells_.size(), comm_->rank());

    globalIds_.resize(globalCells_.size());
    std::iota(globalIds_.begin(), globalIds_.end(), 0);

    bBox_ = BoundingBox(nodes_.begin(), nodes_.end());
}

void FiniteVolumeGrid2D::initPatches(const std::unordered_map<std::string, std::vector<Label>> &patches)
{
    patches_.clear();

    for (const auto &entry: patches)
        createPatchByNodes(entry.first, entry.second);
}

void FiniteVolumeGrid2D::initCommBuffers(const std::vector<Label> &ownership, const std::vector<Label> &globalIds)
{
    if (ownership.size() != cells_.size() || globalIds.size() != cells_.size())
        throw Exception("FiniteVolumeGrid2D",
                        "initCommBuffers",
                        "ownership and global id vector size must match the number of cells.");

    cellOwnership_ = ownership;
    globalIds_ = globalIds;

    sendCellGroups_ = std::vector<CellGroup>(comm_->nProcs());
    bufferCellGroups_ = std::vector<CellGroup>(comm_->nProcs());

    localCells_.clear();
    localCells_.add(cells_.begin(), cells_.end());

    for (const Cell &cell: cells_)
        if (cellOwnership_[cell.id()] != comm_->rank())
        {
            bufferCellGroups_[cellOwnership_[cell.id()]].add(cell);
            localCells_.remove(cell);
        }

    std::vector<std::vector<Label>> recvOrders(comm_->nProcs());

    for (int proc = 0; proc < comm_->nProcs(); ++proc)
    {
        if (proc == comm_->rank())
            continue;

        for (const Cell &cell: bufferCellGroups_[proc])
            recvOrders[proc].push_back(globalIds_[cell.id()]);

        comm_->isend(proc, recvOrders[proc], comm_->rank());
    }

    std::unordered_map<Label, Label> globalToLocalIdMap;

    for (Label id = 0; id < cells_.size(); ++id)
        globalToLocalIdMap[globalIds[id]] = id;

    for (int proc = 0; proc < comm_->nProcs(); ++proc)
    {
        if (proc == comm_->rank())
            continue;

        std::vector<Label> sendOrder(comm_->probeSize<Label>(proc, proc));
        comm_->recv(proc, sendOrder, proc);

        for (Label id: sendOrder)
            sendCellGroups_[proc].add(cells_[globalToLocalIdMap[id]]);
    }

    comm_->waitAll();
}

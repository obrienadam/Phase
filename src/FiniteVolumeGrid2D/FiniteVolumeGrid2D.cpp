#include <numeric>

#include <cgnslib.h>
#include <metis.h>

#include "FiniteVolumeGrid2D.h"

FiniteVolumeGrid2D::FiniteVolumeGrid2D()
        :
        interiorFaces_("InteriorFaces"),
        boundaryFaces_("BoundaryFaces")
{
    auto registry = std::make_shared<CellZone::ZoneRegistry>();
    localActiveCells_ = CellZone("LocalActiveCells", registry);
    localInactiveCells_ = CellZone("GlobalActiveCells", registry);

    registry = std::make_shared<CellZone::ZoneRegistry>();
    globalActiveCells_ = CellZone("GlobalActiveCells", registry);
    globalInactiveCells_ = CellZone("GlobalInactiveCells", registry);

    //- Public registries
    cellZoneRegistry_ = std::make_shared<CellZone::ZoneRegistry>();
    patchRegistry_ = std::make_shared<Patch::PatchRegistry>();

    comm_ = std::make_shared<Communicator>();
}

FiniteVolumeGrid2D::FiniteVolumeGrid2D(const std::vector<Point2D> &nodes,
                                       const std::vector<Label> &cellInds,
                                       const std::vector<Label> &cells)
        :
        FiniteVolumeGrid2D()
{
    init(nodes, cellInds, cells);
}

void FiniteVolumeGrid2D::init(const std::vector<Point2D> &nodes,
                              const std::vector<Label> &cellInds,
                              const std::vector<Label> &cells)
{
    reset();

    for (const Point2D &node: nodes)
        addNode(node);

    for (int i = 0; i < cellInds.size() - 1; ++i)
    {
        int start = cellInds[i];
        int end = cellInds[i + 1];

        std::vector<Label> ids;
        for (int j = start; j < end; ++j)
            ids.push_back(cells[j]);

        createCell(ids);
    }

    initConnectivity();
    computeBoundingBox();
}

void FiniteVolumeGrid2D::reset()
{
    //- Node related data
    nodes_.clear();
    interiorNodes_.clear();
    boundaryNodes_.clear();

    //- Cell related data
    cells_.clear();
    nActiveCellsGlobal_ = 0;

    //- Local cell zones
    localActiveCells_.clear();
    localInactiveCells_.clear();

    //- Cells on this proc that are globally active/inactive
    globalActiveCells_.clear();
    globalInactiveCells_.clear();

    //- User defined cell groups/zones
    cellGroups_.clear();
    cellZones_.clear();

    //- Communication zones
    sendCellGroups_.clear(); // shared pointers are used so that zones can be moveable!
    bufferCellZones_.clear();

    //- Face related data
    faces_.clear();
    faceDirectory_.clear(); // A directory that can find a face given the two node ids

    //- Interior and boundary face data structures
    interiorFaces_.clear();
    boundaryFaces_.clear();

    //- User defined face groups and patches
    faceGroups_.clear();
    patches_.clear();

    bBox_ = BoundingBox(Point2D(0., 0.), Point2D(0., 0.));
    nodeGroup_.clear();
}

//- size info
std::string FiniteVolumeGrid2D::gridInfo() const
{
    using namespace std;

    return "Finite volume grid info:\n"
                   "------------------------\n"
                   "Number of nodes: " + to_string(nodes_.size()) + "\n" +
           "Number of cells total: " + to_string(cells_.size()) + "\n" +
           "Number of interior faces: " + to_string(interiorFaces_.size()) + "\n" +
           "Number of boundary faces: " + to_string(boundaryFaces_.size()) + "\n" +
           "Number of faces total: " + to_string(faces_.size()) + "\n" +
           "Bounding box: " + bBox_.toString() + "\n";
}

//- Create grid entities
Label FiniteVolumeGrid2D::createCell(const std::vector<Label> &nodeIds)
{
    using namespace std;

    cells_.push_back(Cell(nodeIds, *this));
    Cell &newCell = cells_.back();

    for (Label id: nodeIds)
        nodes_[id].addCell(newCell);

    for (Label i = 0, end = nodeIds.size(); i < end; ++i) {
        Label n1 = nodeIds[i], n2 = nodeIds[(i + 1) % end];

        auto key = n1 < n2 ? make_pair(n1, n2) : make_pair(n2, n1);
        auto it = faceDirectory_.find(key);

        if (it == faceDirectory_.end()) // face doesn't exist, so create it
        {
            faces_.push_back(Face(n1, n2, *this, Face::BOUNDARY));

            Face &face = faces_.back();

            faceDirectory_[key] = face.id();

            face.addCell(newCell);
        } else // face already exists, but is now an interior face
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

std::vector<Point2D> FiniteVolumeGrid2D::coords() const
{
    std::vector<Point2D> coords;
    std::transform(nodes_.begin(), nodes_.end(),
                   std::back_inserter(coords),
                   [](const Node &node) { return Point2D(node); });

    return coords;
}

std::vector<Scalar> FiniteVolumeGrid2D::xCoords() const
{
    std::vector<Scalar> xCoords;

    std::transform(nodes_.begin(), nodes_.end(),
                   std::back_inserter(xCoords),
                   [](const Node &node) { return node.x; });

    return xCoords;
}

std::vector<Scalar> FiniteVolumeGrid2D::yCoords() const
{
    std::vector<Scalar> yCoords;

    std::transform(nodes_.begin(), nodes_.end(),
                   std::back_inserter(yCoords),
                   [](const Node &node) { return node.y; });

    return yCoords;
}

std::vector<Ref<Node> > FiniteVolumeGrid2D::getNodes(const std::vector<Label> &ids)
{
    std::vector<Ref<Node>> nodes;
    nodes.reserve(ids.size());

    for (Label id: ids)
        nodes.push_back(nodes_[id]);

    return nodes;
}

std::vector<Ref<const Node> > FiniteVolumeGrid2D::getNodes(const std::vector<Label> &ids) const
{
    std::vector<Ref<const Node>> nodes;
    nodes.reserve(ids.size());

    for (Label id: ids)
        nodes.push_back(nodes_[id]);

    return nodes;
}

void FiniteVolumeGrid2D::assignNodeIds()
{
    Label id = 0;
    for (Node &node: nodes_)
        node.setId(id++);
}

std::vector<int> FiniteVolumeGrid2D::elementList() const
{
    std::vector<int> elems;
    elems.reserve(5 * cells_.size());

    for (const Cell &cell: cells_)
    {
        elems.push_back(cell.nodes().size() == 3 ? CGNS_ENUMV(TRI_3) : CGNS_ENUMV(QUAD_4));

        for (const Node &node: cell.nodes())
            elems.push_back(node.id() + 1);
    }

    return elems;
}

void FiniteVolumeGrid2D::setCellActive(const Cell &cell)
{
    localActiveCells_.add(cell);
    globalActiveCells_.add(cell);
}

void FiniteVolumeGrid2D::setCellInactive(const Cell& cell)
{
    localInactiveCells_.add(cell);
    globalInactiveCells_.add(cell);
}

CellGroup FiniteVolumeGrid2D::globalCellGroup(const CellGroup &localGroup) const
{
    std::vector<int> isInGlobalGroup(cells_.size(), 0);
    CellGroup globalGroup(localGroup.name());

    for(const Cell& cell: localGroup)
    {
        isInGlobalGroup[cell.id()] = 1;
        globalGroup.add(cell);
    }

    sendMessages(isInGlobalGroup);

    for(const CellZone& bufferZone: bufferZones())
        for(const Cell& cell: bufferZone)
            if(isInGlobalGroup[cell.id()])
                globalGroup.add(cell);

    return globalGroup;
}

CellGroup &FiniteVolumeGrid2D::createCellGroup(const std::string& name)
{
    return *(cellGroups_[name] = std::make_shared<CellGroup>(name));
}

CellZone &FiniteVolumeGrid2D::createCellZone(const std::string& name, std::shared_ptr<CellZone::ZoneRegistry> registry)
{
    return *(cellZones_[name] = std::make_shared<CellZone>(name, registry ? registry: cellZoneRegistry_));
}

//- Face related methods

bool FiniteVolumeGrid2D::faceExists(Label n1, Label n2) const
{
    using namespace std;

    auto it = faceDirectory_.find(
            n1 < n2 ? make_pair(n1, n2) : make_pair(n2, n1)
    );

    return it != faceDirectory_.end();
}

Label FiniteVolumeGrid2D::findFace(Label n1, Label n2) const
{
    using namespace std;

    auto it = faceDirectory_.find(
            n1 < n2 ? make_pair(n1, n2) : make_pair(n2, n1)
    );

    if (it == faceDirectory_.end())
        throw Exception("FiniteVolumeGrid2D", "findFace",
                        "no face found between n1 = " + to_string(n1) + ", n2 = " + to_string(n2) + ".");

    return it->second;
}

void FiniteVolumeGrid2D::assignFaceIds()
{
    Label id = 0;
    for (Face &face: faces_)
        face.setId(id++);
}

FaceGroup &FiniteVolumeGrid2D::createFaceGroup(const std::string &name, const std::vector<Label> &ids)
{
    FaceGroup& faceGroup = faceGroups_[name] = FaceGroup(name);

    for(Label id: ids)
        faceGroup.add(faces_[id]);

    return faceGroup;
}

//- Patch related methods
Patch &FiniteVolumeGrid2D::createPatch(const std::string &name, const std::vector<Label> &faces)
{
    Patch& patch = patches_.insert(std::make_pair(name, Patch(name, patches_.size(), patchRegistry_))).first->second;
    patch.clear();

    for (Label id: faces)
        patch.add(faces_[id]);

    return patch;
}

Patch &FiniteVolumeGrid2D::createPatchByNodes(const std::string &name, const std::vector<Label> &nodes)
{
    Patch& patch = patches_.insert(std::make_pair(name, Patch(name, patches_.size(), patchRegistry_))).first->second;
    patch.clear();

    for (Label lid = 0, rid = 1; rid < nodes.size(); lid += 2, rid += 2)
        patch.add(faces_[findFace(nodes[lid], nodes[rid])]);

    return patch;
}

const std::vector<Ref<const Patch> > FiniteVolumeGrid2D::patches() const
{
    std::vector<Ref<const Patch>> patches;

    for(const auto& entry: patches_)
        patches.push_back(std::cref(entry.second));

    return patches;
}

const Node &FiniteVolumeGrid2D::findNearestNode(const Point2D &pt) const
{
    return nodeGroup_.nearestItems(pt, 1)[0];
}

std::vector<Ref<const Node>> FiniteVolumeGrid2D::findNearestNodes(const Point2D& pt, int nNodes) const
{
    return nodeGroup_.nearestItems(pt, nNodes);
}

std::pair<std::vector<int>, std::vector<int> > FiniteVolumeGrid2D::nodeElementConnectivity() const
{
    using namespace std;
    pair<vector<int>, vector<int>> connectivity;
    connectivity.first.push_back(0);

    for (const Cell &cell: cells_)
    {
        connectivity.first.push_back(connectivity.first.back() + cell.nodes().size());

        for (const Node &node: cell.nodes())
            connectivity.second.push_back(node.id());
    }

    return connectivity;
}

void FiniteVolumeGrid2D::partition(const Input &input, std::shared_ptr<Communicator> comm)
{
    using namespace std;

    comm_ = comm;
    if (comm_->nProcs() == 1) // no need to perform a partition
        return;

    comm_->printf("Partitioning grid into %d partitions...\n", comm_->nProcs());
    vector<idx_t> cellPartition(nCells());

    if (comm_->isMainProc()) // partition is performed on main proc
    {
        idx_t nPartitions = comm_->nProcs();
        idx_t nElems = nCells();
        pair<vector<int>, vector<int>> mesh = nodeElementConnectivity();
        idx_t nNodes = this->nNodes();
        idx_t nCommon = 2; //- face connectivity weighting only
        idx_t objVal;
        vector<idx_t> nodePartition(this->nNodes());

        int status = METIS_PartMeshDual(&nElems, &nNodes,
                                        mesh.first.data(), mesh.second.data(),
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
    auto addCellToThisProc = [this, &cellPartition](const Cell& cell, Scalar r = 0.)->bool {
        if(cellPartition[cell.id()] == comm_->rank())
            return true;

        for(const InteriorLink& nb: cell.neighbours())
            if(cellPartition[nb.cell().id()] == comm_->rank())
                return true;

        for(const CellLink& dg: cell.diagonals())
            if(cellPartition[dg.cell().id()] == comm_->rank())
                return true;

        for(const Cell& kCell: globalActiveCells_.itemsWithin(Circle(cell.centroid(), r)))
            if(cellPartition[kCell.id()] == comm_->rank())
                return true;

        return false;
    };

    //- Construct the crs representation of the local grid
    comm_->printf("Computing the local cell domains...\n");
    vector<Point2D> nodes;
    vector<Label> cellInds(1, 0), cellNodeIds, cellProc;
    unordered_map<Label, Label> cellGlobalToLocalIdMap, cellLocalToGlobalIdMap;
    vector<int> localNodeId(nodes_.size(), -1);
    Scalar r = input.caseInput().get<Scalar>("Grid.minBufferWidth", 0.); //- May be important for algorithms requiring spatial searches

    for(const Cell& cell: cells_)
    {
        if(addCellToThisProc(cell, r))
        {
            cellInds.push_back(cellInds.back() + cell.nodes().size());
            cellProc.push_back(cellPartition[cell.id()]);

            cellGlobalToLocalIdMap[cell.id()] = cellInds.size() - 2;
            cellLocalToGlobalIdMap[cellInds.size() - 2] = cell.id();

            for(const Node& node: cell.nodes())
            {
                if(localNodeId[node.id()] == -1)
                {
                    localNodeId[node.id()] = nodes.size();
                    nodes.push_back(node);
                }

                cellNodeIds.push_back(localNodeId[node.id()]);
            }
        }
    }

    //- Boundary patches
    comm_->printf("Computing the local boundary patches...\n");
    std::unordered_map<std::string, std::vector<Label>> localPatches;

    for(const Patch& patch: patches())
    {
        vector<Label> nodeIds;

        for (const Face &face: patch)
        {
            int lid = localNodeId[face.lNode().id()];
            int rid = localNodeId[face.rNode().id()];

            if (lid > -1 && rid > -1)
            {
                nodeIds.push_back((Label)lid);
                nodeIds.push_back((Label)rid);
            }
        }

        if(!nodeIds.empty())
            localPatches[patch.name()] = nodeIds;
    }

    //- Now re-initialize local domains
    comm_->printf("Initializing local domains...\n");
    init(nodes, cellInds, cellNodeIds);
    for(const auto& patch: localPatches)
        createPatchByNodes(patch.first, patch.second);

    comm_->printf("Finished initializing local domains.\n");

    //- Interprocess communication zones
    comm_->printf("Initializing interprocess communication buffers...\n");
    sendCellGroups_.resize(comm_->nProcs());
    bufferCellZones_.resize(comm_->nProcs());

    for(int proc = 0; proc < comm_->nProcs(); ++proc)
    {
        sendCellGroups_[proc] = CellGroup("Proc" + std::to_string(proc));
        bufferCellZones_[proc] = CellZone("Proc" + std::to_string(proc), localActiveCells_.registry());
    }

    //- Identify buffer regions
    for(const Cell& cell: cells_)
        if(cellProc[cell.id()] != comm_->rank())
            bufferCellZones_[cellProc[cell.id()]].add(cell);

    //- Initialize all buffers
    std::vector<std::vector<unsigned long>> recvOrders(comm_->nProcs());

    for(int proc = 0; proc < comm_->nProcs(); ++proc)
    {
        std::transform(bufferCellZones_[proc].begin(), bufferCellZones_[proc].end(),
                       std::back_inserter(recvOrders[proc]), [&cellLocalToGlobalIdMap](const Cell &cell) {
                    return cellLocalToGlobalIdMap[cell.id()];
                });

        comm_->isend(proc, recvOrders[proc], proc);
    }

    for(int proc = 0; proc < comm_->nProcs(); ++proc)
    {
        std::vector<unsigned long> sendOrder(comm_->probeSize<unsigned long>(proc, comm_->rank()));
        comm_->recv(proc, sendOrder, comm_->rank());

        for(Label gid: sendOrder)
            sendCellGroups_[proc].add(cells_[cellGlobalToLocalIdMap[gid]]);
    }

    comm_->waitAll();
    computeGlobalOrdering();
}

void FiniteVolumeGrid2D::computeGlobalOrdering()
{
    /* Note: global orderings are guaranteed to be contigous on each process.
     * This is a requirement for the majority of sparse linear solvers. */

    std::vector<Size> nLocalCells = comm_->allGather(nLocalActiveCells());
    Index globalIndexStart = std::accumulate(nLocalCells.begin(), nLocalCells.begin() + comm_->rank(), 0);
    nActiveCellsGlobal_ = std::accumulate(nLocalCells.begin(), nLocalCells.end(), 0);

    std::vector<Index> globalIndices[] = {
            std::vector<Index>(nCells(), -1),
            std::vector<Index>(nCells(), -1),
            std::vector<Index>(nCells(), -1),
    };

    Index localIndex = 0;
    for (const Cell &cell: localActiveCells_)
    {
        cells_[cell.id()].setNumIndices(4);
        cells_[cell.id()].index(0) = localIndex;
        cells_[cell.id()].index(1) = globalIndexStart + localIndex;
        cells_[cell.id()].index(2) = 2 * globalIndexStart + localIndex;
        cells_[cell.id()].index(3) = 2 * globalIndexStart + localIndex + nLocalCells[comm_->rank()];
        ++localIndex;

        globalIndices[0][cell.id()] = cell.index(1);
        globalIndices[1][cell.id()] = cell.index(2);
        globalIndices[2][cell.id()] = cell.index(3);
    }

    //- Must communicate new global indices to neighbours
    sendMessages(globalIndices[0]);
    sendMessages(globalIndices[1]);
    sendMessages(globalIndices[2]);

    //- Set global ids for the buffer zones
    for (CellZone &bufferZone: bufferCellZones_)
        for (const Cell &cell: bufferZone)
        {
            cells_[cell.id()].setNumIndices(4);
            cells_[cell.id()].index(0) = -1;
            cells_[cell.id()].index(1) = globalIndices[0][cell.id()];
            cells_[cell.id()].index(2) = globalIndices[1][cell.id()];
            cells_[cell.id()].index(3) = globalIndices[2][cell.id()];

            if(globalIndices[0][cell.id()] != -1)
                globalActiveCells_.add(cell);
            else
                globalInactiveCells_.add(cell);
        }

    comm_->printf("Num local cells main proc = %d\nNum global cells = %d\n",
                  nLocalCells[comm_->rank()],
                  nActiveCellsGlobal_);
}

void FiniteVolumeGrid2D::computeParMetisGlobalOrdering()
{

}

//- Protected methods

void FiniteVolumeGrid2D::initNodes()
{
    nodeGroup_.clear();
    nodeGroup_.add(nodes_.begin(), nodes_.end());
}

void FiniteVolumeGrid2D::initCells()
{
    for (const Face &face: faces_)
    {
        if (face.isBoundary())
        {
            Cell &cell = cells_[face.lCell().id()];
            cell.addBoundaryLink(face);

            boundaryFaces_.add(face);
            boundaryNodes_.add(face.lNode());
            boundaryNodes_.add(face.rNode());
        }
        else
        {
            Cell &lCell = cells_[face.lCell().id()];
            Cell &rCell = cells_[face.rCell().id()];

            lCell.addInteriorLink(face, rCell);
            rCell.addInteriorLink(face, lCell);

            interiorFaces_.add(face);
        }
    }

    for(const Face& face: interiorFaces_)
    {
        if (!boundaryNodes_.isInGroup(face.lNode()))
            interiorNodes_.add(face.lNode());
        else if (!boundaryNodes_.isInGroup(face.rNode()))
            interiorNodes_.add(face.rNode());
    }

    //- Initialize diagonal links
    for (Cell &cell: cells_)
    {
        for (const Node &node: cell.nodes())
            for (const Cell &kCell: node.cells())
            {
                if (&cell == &kCell)
                    continue;
                else if (!cellsShareFace(cell, kCell))
                    cell.addDiagonalLink(kCell);
            }
    }

    setCellsActive(cells_.begin(), cells_.end());
}

void FiniteVolumeGrid2D::initConnectivity()
{
    initNodes();
    initCells();
}

void FiniteVolumeGrid2D::computeBoundingBox()
{
    bBox_ = BoundingBox(nodes_.data(), nodes_.size());
}

#include <cgnslib.h>
#include <metis.h>

#include "FiniteVolumeGrid2D.h"
#include "Exception.h"

FiniteVolumeGrid2D::FiniteVolumeGrid2D(Size nNodes, Size nCells, Size nFaces)
        :
        localActiveCells_("LocalActiveCells", *this),
        inactiveCells_("InactiveCells", *this),
        globalActiveCells_("GlobalActiveCells", *this)
{
    nodes_.reserve(nNodes);
    cells_.reserve(nCells);
    faces_.reserve(nFaces);
}

void FiniteVolumeGrid2D::init(const std::vector<Point2D> &nodes, const std::vector<Label> &cellInds,
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
    localActiveCells_.clear();
    globalActiveCells_.clear();
    inactiveCells_.clear();
    cellGroups_.clear();
    cellZones_.clear();

    nodes_.clear();
    nodeSearch_.clear();
    cells_.clear();
    faces_.clear();
    neighbouringProcs_.clear();
    sendCellGroups_.clear();
    bufferCellZones_.clear();
    faceDirectory_.clear();
    interiorFaces_.clear();
    boundaryFaces_.clear();
    patches_.clear();
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

void FiniteVolumeGrid2D::setCellActive(Cell &cell)
{
    inactiveCells_.remove(cell);
    localActiveCells_.push_back(cell);
    globalActiveCells_.push_back(cell);
}

void FiniteVolumeGrid2D::setCellsActive(const std::vector<Label> &ids)
{
    for (Cell &cell: getCells(ids))
        setCellActive(cell);
}

void FiniteVolumeGrid2D::setCellsActive(CellGroup &cellGroup)
{
    for (Cell &cell: cellGroup)
        setCellActive(cell);
}

void FiniteVolumeGrid2D::setCellsInactive(const std::vector<Label> &ids)
{
    for (Cell &cell: getCells(ids))
    {
        localActiveCells_.remove(cell);
        globalActiveCells_.remove(cell);
        inactiveCells_.push_back(cell);
        cell.clearIndices();
    }
}

void FiniteVolumeGrid2D::setCellsInactive(CellGroup &cellGroup)
{
    for (Cell &cell: cellGroup)
    {
        localActiveCells_.remove(cell);
        globalActiveCells_.remove(cell);
        inactiveCells_.push_back(cell);
        cell.clearIndices();
    }
}

void FiniteVolumeGrid2D::setCellsLocallyInactive(const std::vector<Label> &ids)
{
    for (const Cell &cell: getCells(ids))
        localActiveCells_.remove(cell);
}

CellGroup &FiniteVolumeGrid2D::createCellGroup(const std::string &name, const std::vector<Label> &ids)
{
    cellGroups_.insert(std::make_pair(name, CellGroup(name, *this)));
    CellGroup &cellGroup = cellGroups_.find(name)->second;

    for (Cell &cell: getCells(ids))
        cellGroup.push_back(cell);

    return cellGroup;
}

CellZone &FiniteVolumeGrid2D::createCellZone(const std::string &name, const std::vector<Label> &ids)
{
    cellZones_.insert(std::make_pair(name, CellZone(name, *this)));
    CellZone &cellZone = cellZones_.find(name)->second;

    for (Cell &cell: getCells(ids))
        cellZone.moveToGroup(cell);

    return cellZone;
}

void FiniteVolumeGrid2D::setAllCellsActive()
{
    inactiveCells_.clear();
    localActiveCells_.clear();
    globalActiveCells_.clear();

    for (Cell &cell: cells_)
    {
        localActiveCells_.push_back(cell);
        globalActiveCells_.push_back(cell);
    }
}

//- Cell related methods

std::vector<Ref<Cell> > FiniteVolumeGrid2D::getCells(const std::vector<Label> &ids)
{
    std::vector<Ref<Cell> > cells;
    cells.reserve(ids.size());

    std::transform(ids.begin(),
                   ids.end(),
                   std::back_inserter(cells),
                   [this](Label id) { return std::ref(this->cells_[id]); });

    return cells;
}

std::vector<Ref<const Cell> > FiniteVolumeGrid2D::getCells(const std::vector<Label> &ids) const
{
    std::vector<Ref<const Cell> > cells;
    cells.reserve(ids.size());

    std::transform(ids.begin(),
                   ids.end(),
                   std::back_inserter(cells),
                   [this](Label id) { return std::cref(this->cells_[id]); });

    return cells;
}

void FiniteVolumeGrid2D::assignCellIds()
{
    Label id = 0;
    for (Cell &cell: cells_)
        cell.setId(id++);
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

//- Patch related methods
void FiniteVolumeGrid2D::applyPatch(const std::string &patchName, const std::vector<Label> &faces)
{
    auto insert = patches_.insert(std::make_pair(patchName, Patch(patchName, patches_.size())));

    if (!insert.second)
        throw Exception("FiniteVolumeGrid2D", "applyPatch", "patch already exists.");

    Patch &patch = (insert.first)->second;

    for (Label id: faces)
        patch.addFace(faces_[id]);
}

void FiniteVolumeGrid2D::applyPatchByNodes(const std::string &patchName, const std::vector<Label> &nodes)
{
    auto insert = patches_.insert(std::make_pair(patchName, Patch(patchName, patches_.size())));

    if (!insert.second)
        throw Exception("FiniteVolumeGrid2D", "applyPatchByNodes", "patch already exists.");

    Patch &patch = (insert.first)->second;

    for (Label lid = 0, rid = 1; rid < nodes.size(); lid += 2, rid += 2)
        patch.addFace(faces_[findFace(nodes[lid], nodes[rid])]);
}

const Node &FiniteVolumeGrid2D::findNearestNode(const Point2D &pt) const
{
    return nodes_[nodeSearch_.kNearestNeighbourSearch(pt, 1)[0]];
}

std::vector<Ref<Node> > FiniteVolumeGrid2D::nodesInShape(const Shape2D &shape)
{
    return getNodes(nodeSearch_.rangeSearch(shape));
}

std::vector<Ref<const Node> > FiniteVolumeGrid2D::nodesInShape(const Shape2D &shape) const
{
    return getNodes(nodeSearch_.rangeSearch(shape));
}

std::vector<std::vector<Ref<const Cell>>> FiniteVolumeGrid2D::constructSmoothingKernels(Scalar width) const
{
    std::vector<std::vector<Ref<const Cell>>> kernels(nCells());
    const CellGroup &group = globalActiveCells();

    for (const Cell &cell: group)
    {
        Circle base = Circle(cell.centroid(), width);
        std::vector<Ref<const Cell>> kernel = group.cellCentersWithin(base);

        //- must  find intersecting boundary edges

        kernels[cell.id()] = kernel;
    }

    return kernels;
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

void FiniteVolumeGrid2D::partition(const Input &input, const Communicator &comm)
{
    using namespace std;
    comm.printf("Partitioning...\n");
    if (comm.nProcs() == 1) // no need to perform a partition
        return;

    vector<idx_t> cellPartition(nCells());

    if (comm.isMainProc())
    {
        idx_t nPartitions = comm.nProcs();
        idx_t nElems = nCells();
        pair<vector<int>, vector<int>> mesh = nodeElementConnectivity();
        idx_t nNodes = this->nNodes();
        idx_t nCommon = 2;
        idx_t objVal;
        vector<idx_t> nodePartition(this->nNodes());

        comm.printf("Partitioning mesh into %d subdomains...\n", nPartitions);
        int status = METIS_PartMeshDual(&nElems, &nNodes,
                                        mesh.first.data(), mesh.second.data(),
                                        NULL, NULL,
                                        &nCommon, &nPartitions,
                                        NULL, NULL, &objVal,
                                        cellPartition.data(), nodePartition.data());
        if (status == METIS_OK)
            comm.printf("Sucessfully partitioned mesh.\n");
        else
            throw Exception("FiniteVolumeGrid2D", "partition", "an error occurred during partitioning.");
    }

    comm.broadcast(comm.mainProcNo(), cellPartition);

    vector<Point2D> localNodes;
    vector<Label> localCellInds(1, 0), localCells, localCellList;
    map<string, vector<Label>> localPatches;
    vector<int> localNodeId(nodes_.size(), -1);
    vector<int> localCellId(cells_.size(), -1);
    Label nextLocalNodeId = 0, nextLocalCellId = 0;

    comm.printf("Initializing local domains...\n");
    auto addCellToPart = [&](const Cell &cell) -> bool // helper lambda, returns whether or not a cell was added
    {
        if (localCellId[cell.id()] != -1)
            return false;

        localCellList.push_back(cell.id());
        localCellInds.push_back(localCellInds.back() + cell.nodes().size());

        for (const Node &node: cell.nodes())
        {
            if (localNodeId[node.id()] == -1)
            {
                localNodes.push_back(nodes_[node.id()]);
                localNodeId[node.id()] = nextLocalNodeId++;
            }

            localCells.push_back(localNodeId[node.id()]);
        }

        //- check for patches
        for (const BoundaryLink &bd: cell.boundaries())
        {
            localPatches[bd.face().patch().name].push_back(bd.face().lNode().id());
            localPatches[bd.face().patch().name].push_back(bd.face().rNode().id());
        }

        localCellId[cell.id()] = nextLocalCellId++;

        return true;
    };

    for (const Cell &cell: cells_)
    {
        if (comm.rank() == cellPartition[cell.id()])
            addCellToPart(cell);
    }

    comm.printf("Finished initializing local domains.\n"
                        "Constructing intraprocess communication buffers...\n");

    vector<int> neighboursProc(comm.nProcs(), -1);
    vector<int> nbProcs;
    vector<vector<Label>> sendOrder, recvOrder;

    int bufferWidth = input.caseInput().get<int>("Grid.partitionBufferWidth", 1); //- Must be set carefully!

    for (int i = 0; i < bufferWidth; ++i) // This is where buffer width can be set
        for (const Cell &cell: getCells(localCellList))
        {
            for (const InteriorLink &nb: cell.neighbours())
            {
                int nbProc = cellPartition[nb.cell().id()];

                if (comm.rank() == nbProc)
                    continue;
                else if (neighboursProc[nbProc] == -1)
                {
                    neighboursProc[nbProc] = nbProcs.size();
                    nbProcs.push_back(nbProc);
                    recvOrder.push_back(vector<Label>());
                }

                if (addCellToPart(nb.cell()))
                    recvOrder[neighboursProc[nbProc]].push_back(nb.cell().id());
            }

            for (const DiagonalCellLink &dg: cell.diagonals())
            {
                int nbProc = cellPartition[dg.cell().id()];

                if (comm.rank() == nbProc)
                    continue;
                else if (neighboursProc[nbProc] == -1)
                {
                    neighboursProc[nbProc] = nbProcs.size();
                    nbProcs.push_back(nbProc);
                    recvOrder.push_back(vector<Label>());
                }

                if (addCellToPart(dg.cell()))
                    recvOrder[neighboursProc[nbProc]].push_back(dg.cell().id());
            }
        }

    //- Send and receive the buffer orders
    vector<vector<int>> sendSizes(nbProcs.size(), vector<int>(1));
    for (int i = 0; i < nbProcs.size(); ++i)
        comm.irecv(nbProcs[i], sendSizes[i]);

    for (int i = 0; i < nbProcs.size(); ++i)
        comm.ssend(nbProcs[i], vector<int>(1, recvOrder[i].size()));

    sendOrder.resize(nbProcs.size());
    comm.waitAll();

    for (int i = 0; i < nbProcs.size(); ++i)
    {
        sendOrder[i].resize(sendSizes[i][0]);
        comm.irecv(nbProcs[i], sendOrder[i]);
    }

    for (int i = 0; i < nbProcs.size(); ++i)
        comm.ssend(nbProcs[i], recvOrder[i]);

    comm.waitAll();

    comm.printf("Finished constructing intraprocess communication buffers.\n"
                        "Re-initializing local grid entities...\n");

    init(localNodes, localCellInds, localCells);

    //- Map the buffer orders to local ids
    for (int i = 0; i < nbProcs.size(); ++i)
    {
        for (Label &id: sendOrder[i])
            id = localCellId[id];

        for (Label &id: recvOrder[i])
            id = localCellId[id];

        addNeighbouringProc(nbProcs[i], sendOrder[i], recvOrder[i]);
    }

    comm.printf("Finished re-initializing local grid entities.\n"
                        "Initialing local patches...\n");

    //- Map boundary patches to local ids
    for (auto &entry: localPatches)
    {
        for (Label &id: entry.second)
            id = localNodeId[id];

        applyPatchByNodes(entry.first, entry.second);
    }

    vector<Label> partitionPatch;
    for (const Face &face: boundaryFaces())
        if (!face.belongsToPatch())
            partitionPatch.push_back(face.id());

    applyPatch("LocalCellBoundary", partitionPatch);
    comm.printf("Finished initializing local patches.\n");
    comm.printf("Finished partitioning grid.\n");
}

void FiniteVolumeGrid2D::addNeighbouringProc(int procNo,
                                             const std::vector<Label> &sendOrder,
                                             const std::vector<Label> &recvOrder)
{
    neighbouringProcs_.push_back(procNo);

    CellGroup sendCellGroup("SendCellGroupProc" + std::to_string(procNo), *this);

    for (Label id: sendOrder)
        sendCellGroup.push_back(cells_[id]);

    CellZone bufferCellZone("RecvCellZoneProc" + std::to_string(procNo), *this);

    for (Label id: recvOrder)
        bufferCellZone.push_back(cells_[id]);

    sendCellGroups_.push_back(sendCellGroup);
    bufferCellZones_.push_back(bufferCellZone);

    setCellsLocallyInactive(recvOrder);
}

void FiniteVolumeGrid2D::computeGlobalOrdering(const Communicator &comm)
{
    /* Note: global orderings are guaranteed to be contigous on each process.
     * This is a requirement for the majority of sparse linear solvers. */

    std::vector<Size> nLocalCells = comm.allGather(nLocalActiveCells());

    Index globalIndexStart = 0;
    for (int proc = 0; proc < comm.rank(); ++proc)
        globalIndexStart += nLocalCells[proc];

    std::vector<Index> globalIndices[] = {
            std::vector<Index>(nCells(), -1),
            std::vector<Index>(nCells(), -1),
            std::vector<Index>(nCells(), -1),
    };

    Index localIndex = 0;
    for (Cell &cell: localActiveCells_)
    {
        cell.setNumIndices(4);
        cell.index(0) = localIndex;
        cell.index(1) = globalIndexStart + localIndex;
        cell.index(2) = 2 * globalIndexStart + localIndex;
        cell.index(3) = 2 * globalIndexStart + localIndex + nLocalCells[comm.rank()];
        ++localIndex;

        globalIndices[0][cell.id()] = cell.index(1);
        globalIndices[1][cell.id()] = cell.index(2);
        globalIndices[2][cell.id()] = cell.index(3);
    }

    //- Must communicate new global indices to neighbours
    sendMessages(comm, globalIndices[0]);
    sendMessages(comm, globalIndices[1]);
    sendMessages(comm, globalIndices[2]);

    //- Set global ids for the buffer zones
    for (const CellZone &bufferZone: bufferCellZones_)
        for (Cell &cell: bufferZone)
        {
            cell.setNumIndices(4);
            cell.index(0) = -1;
            cell.index(1) = globalIndices[0][cell.id()];
            cell.index(2) = globalIndices[1][cell.id()];
            cell.index(3) = globalIndices[2][cell.id()];
        }

    nActiveCellsGlobal_ = comm.sum(nLocalActiveCells());
    comm.printf("Num local cells main proc = %d\nNum global cells = %d\n", nLocalCells[comm.rank()],
                nActiveCellsGlobal_);
}

//- Protected methods

void FiniteVolumeGrid2D::initNodes()
{
    nodeSearch_.constructRTree(coords());
}

void FiniteVolumeGrid2D::initCells()
{
    for (const Face &face: faces_)
    {
        if (face.isBoundary())
        {
            Cell &cell = cells_[face.lCell().id()];
            cell.addBoundaryLink(face);

            boundaryFaces_.push_back(face);
        }
        else
        {
            Cell &lCell = cells_[face.lCell().id()];
            Cell &rCell = cells_[face.rCell().id()];

            lCell.addInteriorLink(face, rCell);
            rCell.addInteriorLink(face, lCell);

            interiorFaces_.push_back(face);
        }
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

    setAllCellsActive();
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

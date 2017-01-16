#include <cgnslib.h>
#include <metis.h>

#include "FiniteVolumeGrid2D.h"
#include "Exception.h"

FiniteVolumeGrid2D::FiniteVolumeGrid2D(Size nNodes, Size nCells, Size nFaces)
    :
      nodeSearch_(nodes_)
{
    nodes_.reserve(nNodes);
    cells_.reserve(nCells);
    faces_.reserve(nFaces);

    cellGroups_["active"] = CellGroup("active");
    cellZones_["fluid"] = CellZone("fluid");
    cellZones_["inactive"] = CellZone("inactive");
}

void FiniteVolumeGrid2D::init(const std::vector<Point2D> &nodes, const std::vector<Label> &cellInds, const std::vector<Label> &cells)
{
    reset();

    for(const Point2D& node: nodes)
        addNode(node);

    for(int i = 0; i < cellInds.size() - 1; ++i)
    {
        int start = cellInds[i];
        int end = cellInds[i + 1];
        std::vector<Label> ids;

        for(int j = start; j < end; ++j)
            ids.push_back(cells[j]);

        createCell(ids);
    }

    initConnectivity();
    computeBoundingBox();
}

void FiniteVolumeGrid2D::reset()
{
    nodes_.clear();
    cells_.clear();
    faces_.clear();;
    cellGroups_.clear();
    cellZones_.clear();
    neighbouringProcs_.clear();
    procSendOrder_.clear();
    procRecvOrder_.clear();
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
Label FiniteVolumeGrid2D::createCell(const std::vector<Label>& nodeIds)
{
    using namespace std;

    cells_.push_back(Cell(nodeIds, *this));

    Cell &newCell = cells_.back();

    for(Label id: nodeIds)
        nodes_[id].addCell(newCell);

    for(Label i = 0, end = nodeIds.size(); i < end; ++i)
    {
        Label n1 = nodeIds[i], n2 = nodeIds[(i + 1)%end];

        auto key = n1 < n2 ? make_pair(n1, n2) : make_pair(n2, n1);
        auto it = faceDirectory_.find(key);

        if(it == faceDirectory_.end()) // face doesn't exist, so create it
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
                   [](const Node& node){ return Point2D(node); });

    return coords;
}

std::vector<Scalar> FiniteVolumeGrid2D::xCoords() const
{
    std::vector<Scalar> xCoords;

    std::transform(nodes_.begin(), nodes_.end(),
                   std::back_inserter(xCoords),
                   [](const Node& node){ return node.x; });

    return xCoords;
}

std::vector<Scalar> FiniteVolumeGrid2D::yCoords() const
{
    std::vector<Scalar> yCoords;

    std::transform(nodes_.begin(), nodes_.end(),
                   std::back_inserter(yCoords),
                   [](const Node& node){ return node.y; });

    return yCoords;
}

void FiniteVolumeGrid2D::assignNodeIds()
{
    Label id = 0;
    for(Node &node: nodes_)
        node.setId(id++);
}

std::vector<int> FiniteVolumeGrid2D::elementList() const
{
    std::vector<int> elems;
    elems.reserve(5*cells_.size());

    for(const Cell& cell: cells_)
    {
        elems.push_back(cell.nodes().size() == 3 ? TRI_3 : QUAD_4);

        for(const Node& node: cell.nodes())
            elems.push_back(node.id() + 1);
    }

    return elems;
}

CellZone &FiniteVolumeGrid2D::addCellZone(const std::string &name, const std::vector<Label> &ids)
{
    CellZone& newCellZone = cellZones_[name];
    newCellZone.moveToGroup(getCells(ids));
    return newCellZone;
}

//- Cell related methods

CellZone& FiniteVolumeGrid2D::moveCellsToFluidCellGroup(const std::vector<size_t>& ids)
{
    CellZone& fluidCells = cellZones_.find("fluid")->second;

    for(size_t id: ids)
    {
        cells_[id].setActive();
        cells_[id].setFluidCell();
        fluidCells.moveToGroup(cells_[id]);
    }

    constructActiveCellGroup();
    return fluidCells;
}

CellZone& FiniteVolumeGrid2D::moveAllCellsToFluidCellGroup()
{
    for(const Cell &cell: cells_)
    {
        cell.setActive();
        cell.setFluidCell();
    }

    cellZones_.find("fluid")->second.moveAllCellsToThisGroup();

    constructActiveCellGroup();
    return cellZones_.find("fluid")->second;
}

CellZone& FiniteVolumeGrid2D::moveCellsToInactiveCellGroup(const std::vector<size_t>& ids)
{
    CellZone& inactiveCells = cellZones_.find("inactive")->second;

    for(size_t id: ids)
    {
        cells_[id].setInactive();
        cells_[id].setNonFluidCell();
        inactiveCells.moveToGroup(cells_[id]);
    }

    constructActiveCellGroup();
    return inactiveCells;
}

CellZone& FiniteVolumeGrid2D::moveCellsToCellGroup(const std::string& name, const std::vector<size_t>& ids)
{
    CellZone &group = (cellZones_.insert(std::make_pair(name, CellZone(name)))).first->second;

    for(size_t id: ids)
    {
        cells_[id].setActive();
        cells_[id].setNonFluidCell();
        group.moveToGroup(cells_[id]);
    }

    constructActiveCellGroup();
    return group;
}

void FiniteVolumeGrid2D::removeFromActiveCellGroup(const std::vector<Label> &ids)
{
    CellGroup& group = cellGroups_["active"];

    for(const Cell& cell: getCells(ids))
        group.remove(cell);
}

const std::vector<Ref<const Cell> > FiniteVolumeGrid2D::getCells(const std::vector<Label> &ids) const
{
    std::vector<Ref<const Cell> > cells;
    cells.reserve(ids.size());

    std::transform(ids.begin(),
                   ids.end(),
                   std::back_inserter(cells),
                   [this](Label id){ return std::cref(this->cells_[id]); });

    return cells;
}

void FiniteVolumeGrid2D::assignCellIds()
{
    Label id = 0;
    for(Cell &cell: cells_)
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

    if(it == faceDirectory_.end())
        throw Exception("FiniteVolumeGrid2D", "findFace", "no face found between n1 = " + to_string(n1) + ", n2 = " + to_string(n2) + ".");

    return it->second;
}

void FiniteVolumeGrid2D::assignFaceIds()
{
    Label id = 0;
    for(Face &face: faces_)
        face.setId(id++);
}

//- Patch related methods
void FiniteVolumeGrid2D::applyPatch(const std::string &patchName, const std::vector<Label> &faces)
{
    auto insert = patches_.insert(std::make_pair(patchName, Patch(patchName, patches_.size())));

    if(!insert.second)
        throw Exception("FiniteVolumeGrid2D", "applyPatch", "patch already exists.");

    Patch &patch = (insert.first)->second;

    for(Label id: faces)
        patch.addFace(faces_[id]);
}

void FiniteVolumeGrid2D::applyPatchByNodes(const std::string &patchName, const std::vector<Label> &nodes)
{
    auto insert = patches_.insert(std::make_pair(patchName, Patch(patchName, patches_.size())));

    if(!insert.second)
        throw Exception("FiniteVolumeGrid2D", "applyPatchByNodes", "patch already exists.");

    Patch &patch = (insert.first)->second;

    for(Label lid = 0, rid = 1; rid < nodes.size(); lid += 2, rid += 2)
        patch.addFace(faces_[findFace(nodes[lid], nodes[rid])]);
}

const Node& FiniteVolumeGrid2D::findNearestNode(const Point2D& pt) const
{
    return nodeSearch_.kNearestNeighbourSearch(pt, 1)[0];
}

std::vector<std::vector<Ref<const Cell>> > FiniteVolumeGrid2D::constructSmoothingKernels(Scalar width) const
{
    std::vector<std::vector<Ref<const Cell>>> kernels(nCells());
    const CellGroup& group = activeCells();
    BoundingBox gridgeom = boundingBox();

    for(const Cell& cell: group)
    {
        Circle base = Circle(cell.centroid(), width);
        std::vector<Ref<const Cell>> kernel = group.rangeSearch(base);

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

    for(const Cell& cell: cells_)
    {
        connectivity.first.push_back(connectivity.first.back() + cell.nodes().size());

        for(const Node& node: cell.nodes())
            connectivity.second.push_back(node.id());
    }

    return connectivity;
}

void FiniteVolumeGrid2D::partition(const Communicator &comm)
{
    using namespace std;

    if(comm.nProcs() == 1) // no need to perform a partition
        return;

    vector<idx_t> cellPartition(nCells());

    if(comm.isMainProc())
    {
        idx_t nPartitions = comm.nProcs();
        idx_t nElems = nCells();
        pair<vector<int>, vector<int>> mesh = nodeElementConnectivity();
        idx_t nNodes = this->nNodes();
        idx_t nCommon = 2;
        idx_t objVal;
        vector<idx_t> nodePartition(this->nNodes());

        printf("Partitioning mesh into %d subdomains...\n", nPartitions);
        int status = METIS_PartMeshDual(&nElems, &nNodes,
                                        mesh.first.data(), mesh.second.data(),
                                        NULL, NULL,
                                        &nCommon, &nPartitions,
                                        NULL, NULL, &objVal,
                                        cellPartition.data(), nodePartition.data());
        if(status == METIS_OK)
            printf("Sucessfully partitioned mesh.\n");
        else
            throw Exception("finiteVolumeGrid2D", "partition", "an error occurred during partitioning.");
    }

    comm.broadcast(comm.mainProcNo(), cellPartition);

    vector<Point2D> localNodes;
    vector<Label> localCellInds(1, 0), localCells, localCellList;
    map<string, vector<Label>> localPatches;
    vector<int> localNodeId(nodes_.size(), -1);
    vector<int> localCellId(cells_.size(), -1);
    Label nextLocalNodeId = 0, nextLocalCellId = 0;

    auto addCellToPart = [&](const Cell& cell)->bool // helper lambda, returns whether or not a cell was added
    {
        if(localCellId[cell.id()] != -1)
            return false;

        localCellList.push_back(cell.id());
        localCellInds.push_back(localCellInds.back() + cell.nodes().size());

        for(const Node& node: cell.nodes())
        {
            if(localNodeId[node.id()] == -1)
            {
                localNodes.push_back(nodes_[node.id()]);
                localNodeId[node.id()] = nextLocalNodeId++;
            }

            localCells.push_back(localNodeId[node.id()]);
        }

        //- check for patches
        for(const BoundaryLink& bd: cell.boundaries())
        {
            localPatches[bd.face().patch().name].push_back(bd.face().lNode().id());
            localPatches[bd.face().patch().name].push_back(bd.face().rNode().id());
        }

        localCellId[cell.id()] = nextLocalCellId++;

        return true;
    };

    for(const Cell& cell: cells_)
    {
        if(comm.rank() == cellPartition[cell.id()])
            addCellToPart(cell);
    }

    vector<int> neighboursProc(comm.nProcs(), -1);
    vector<int> nbProcs;
    vector<vector<Label>> sendOrder, recvOrder;
for(int i = 0; i < 2; ++i)
    for(const Cell& cell: getCells(localCellList))
    {
        for(const InteriorLink& nb: cell.neighbours())
        {
            int nbProc = cellPartition[nb.cell().id()];

            if(comm.rank() == nbProc)
                continue;
            else if(neighboursProc[nbProc] == -1)
            {
                neighboursProc[nbProc] = nbProcs.size();
                nbProcs.push_back(nbProc);
                recvOrder.push_back(vector<Label>());
            }

            if(addCellToPart(nb.cell()))
               recvOrder[neighboursProc[nbProc]].push_back(nb.cell().id());
        }

        for(const DiagonalCellLink& dg: cell.diagonals())
        {
            int nbProc = cellPartition[dg.cell().id()];

            if(comm.rank() == nbProc)
                continue;
            else if(neighboursProc[nbProc] == -1)
            {
                neighboursProc[nbProc] = nbProcs.size();
                nbProcs.push_back(nbProc);
                recvOrder.push_back(vector<Label>());
            }

            if(addCellToPart(dg.cell()))
                recvOrder[neighboursProc[nbProc]].push_back(dg.cell().id());
        }
    }

    //- Send and receive the buffer orders
    vector<vector<int>> sendSizes(nbProcs.size(), vector<int>(1));
    for(int i = 0; i < nbProcs.size(); ++i)
        comm.irecv(nbProcs[i], sendSizes[i]);

    for(int i = 0; i < nbProcs.size(); ++i)
        comm.send(nbProcs[i], vector<int>(1, recvOrder[i].size()));

    sendOrder.resize(nbProcs.size());
    comm.waitAll();

    for(int i = 0; i < nbProcs.size(); ++i)
    {
        sendOrder[i].resize(sendSizes[i][0]);
        comm.irecv(nbProcs[i], sendOrder[i]);
        comm.send(nbProcs[i], recvOrder[i]);
    }
    comm.waitAll();

    init(localNodes, localCellInds, localCells);

    //- Map the buffer orders to local ids
    for(int i = 0; i < nbProcs.size(); ++i)
    {
        for(Label &id: sendOrder[i])
            id = localCellId[id];

        for(Label &id: recvOrder[i])
            id = localCellId[id];

        addNeighbouringProc(nbProcs[i], sendOrder[i], recvOrder[i]);
    }

    //- Map boundary patches to local ids
    for(auto& entry: localPatches)
    {
        for(Label& id: entry.second)
            id = localNodeId[id];

        applyPatchByNodes(entry.first, entry.second);
    }

    vector<Label> partitionPatch;
    for(const Face& face: boundaryFaces())
        if(!face.belongsToPatch())
            partitionPatch.push_back(face.id());

    applyPatch("partition", partitionPatch);
}

void FiniteVolumeGrid2D::addNeighbouringProc(int procNo,
                                             const std::vector<Label> &sendOrder,
                                             const std::vector<Label> &recvOrder)
{
    neighbouringProcs_.push_back(procNo);
    procSendOrder_.push_back(sendOrder);
    procRecvOrder_.push_back(recvOrder);
    addCellZone("proc" + std::to_string(procNo), recvOrder);
    removeFromActiveCellGroup(recvOrder);
}

std::pair<int, int> FiniteVolumeGrid2D::globalCellRange() const
{
    return std::make_pair(cells_.front().globalId(), cells_.back().globalId() + 1);
}

//- Protected methods

void FiniteVolumeGrid2D::initNodes()
{
    nodeSearch_.constructRTree();
}

void FiniteVolumeGrid2D::initCells()
{
    CellZone& fluidCells = cellZones_["fluid"];

    for(Cell &cell: cells_)
    {
        cell.setActive();
        fluidCells.push_back(cell);
    }

    for(const Face& face: faces_)
    {
        if(face.isBoundary())
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
    for(Cell& cell: cells_)
        for(const Node& node: cell.nodes())
            for(const Cell& kCell: node.cells())
            {
                if(&cell == &kCell)
                    continue;
                else if(!cellsShareFace(cell, kCell))
                    cell.addDiagonalLink(kCell);
            }

    for(Cell& cell: cells_)
    {
        for(const Node& node: cell.nodes())
        {
            for(const Cell& nbCell: node.cells())
            {
                if(cell.id() == nbCell.id())
                    continue;
                else if(!cellsShareFace(cell, nbCell))
                    cell.addDiagonalLink(nbCell);
            }
        }
    }

    constructActiveCellGroup();
}

void FiniteVolumeGrid2D::initConnectivity()
{
    initNodes();
    initCells();
}

void FiniteVolumeGrid2D::initGlobalIds(const Communicator &comm)
{
    std::vector<unsigned long> nCellsPerProc = comm.allGather(cells_.size());
    int lowerId = 0;

    for(int procNo = 0; procNo < comm.rank(); ++procNo)
        lowerId += nCellsPerProc[procNo];

    for(Cell& cell: cells_)
        cell.setGlobalId(lowerId++);
}

void FiniteVolumeGrid2D::constructActiveCellGroup()
{
    Index idx = 0;
    CellGroup& activeGroup = cellGroups_["active"];
    CellZone& inactiveZone = cellZones_["inactive"];

    activeGroup.clear();
    inactiveZone.clear();

    for(Cell &cell: cells_)
    {
        if(cell.isActive())
        {
            activeGroup.push_back(cell);
            cell.setGlobalIndex(idx++);
        }
        else
        {
            inactiveZone.push_back(cell);
            cell.setGlobalIndex(Cell::INACTIVE);
        }
    }
}

void FiniteVolumeGrid2D::computeBoundingBox()
{
    bBox_ = BoundingBox(nodes_.data(), nodes_.size());
}

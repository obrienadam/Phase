#include <cgnslib.h>

#include "FiniteVolumeGrid2D.h"
#include "metis.h"
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

const std::vector<Ref<const Cell> > FiniteVolumeGrid2D::getCells(const std::vector<Label> &ids) const
{
    std::vector<Ref<const Cell> > cells;
    cells.reserve(ids.size());

    for(Label id: ids)
        cells.push_back(std::cref(cells_[id]));

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
        throw Exception("FiniteVolumeGrid2D", "findFace", "face was not found.");

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

const Node& FiniteVolumeGrid2D::findNearestNode(const Point2D& pt) const
{
    return nodeSearch_.kNearestNeighbourSearch(pt, 1)[0];
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

    if(comm.nProcs() == 1)
        return;

    if(comm.mainProc())
    {
        idx_t nPartitions = comm.nProcs();
        idx_t nElems = nCells();
        idx_t nNodes = this->nNodes();
        pair<vector<idx_t>, vector<idx_t>> mesh = nodeElementConnectivity();
        idx_t nCommon = 2;
        idx_t objVal;
        vector<idx_t> cellPartition(nCells());
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

        //- Main partitioning
        vector<vector<Label>> cellLists(comm.nProcs());
        for(Label id = 0; id < cellPartition.size(); ++id)
            cellLists[cellPartition[id]].push_back(id);

        vector<vector<vector<Label>>> sendOrdering(comm.nProcs(), vector<vector<Label>>(comm.nProcs())),
                recvOrdering(comm.nProcs(), vector<vector<Label>>(comm.nProcs()));

        //- Buffer regions
        for(int proc = 0; proc < comm.nProcs(); ++proc)
        {
            vector<bool> isInLocalGroup(nElems, false);

            for(const Cell& cell: getCells(cellLists[proc]))
            {
                for(const InteriorLink& nb: cell.neighbours())
                {
                    int nbProc = cellPartition[nb.cell().id()];

                    if(proc == nbProc)
                        continue;

                    if(!isInLocalGroup[nb.cell().id()])
                    {
                        cellLists[proc].push_back(nb.cell().id());
                        sendOrdering[nbProc][proc].push_back(nb.cell().id());
                        recvOrdering[proc][nbProc].push_back(nb.cell().id());
                        isInLocalGroup[nb.cell().id()] = true;
                    }
                }

                for(const DiagonalCellLink& dg: cell.diagonals())
                {
                    int nbProc = cellPartition[dg.cell().id()];

                    if(proc == nbProc)
                        continue;

                    if(!isInLocalGroup[dg.cell().id()])
                    {
                        cellLists[proc].push_back(dg.cell().id());
                        sendOrdering[nbProc][proc].push_back(dg.cell().id());
                        recvOrdering[proc][nbProc].push_back(dg.cell().id());
                        isInLocalGroup[dg.cell().id()] = true;
                    }
                }
            }
        } // end buffer regions

        //- Figure out neighbouring proces
        vector<vector<int>> neighbouringProcs(comm.nProcs());
        for(int proc = 0; proc < comm.nProcs(); ++proc)
            for(int nbProc = 0; nbProc < comm.nProcs(); ++nbProc)
            {
                if(sendOrdering[proc][nbProc].size() != 0 && recvOrdering[proc][nbProc].size() != 0)
                    neighbouringProcs[proc].push_back(nbProc);
            }

        //- construct local node element lists
        vector<vector<Label>> localCellInds(comm.nProcs()), localCells(comm.nProcs());
        vector<vector<Point2D>> localNodes(comm.nProcs());
        for(int proc = 0; proc < comm.nProcs(); ++proc)
        {
            vector<Label>& cellInds = localCellInds[proc];
            vector<Label>& cells = localCells[proc];
            vector<Point2D>& nodes = localNodes[proc];

            vector<Label> nodeIdMap(nNodes, -1);

            cellInds.push_back(0);

            Label nextLocalNodeId = 0;
            for(const Cell& cell: getCells(cellLists[proc]))
            {
                cellInds.push_back(cellInds.back() + cell.nodes().size());

                for(const Node& node: cell.nodes())
                {
                    if(nodeIdMap[node.id()] == -1)
                    {
                        nodeIdMap[node.id()] = nextLocalNodeId++;
                        nodes.push_back(nodes_[node.id()]);
                    }
                    cells.push_back(nodeIdMap[node.id()]);
                }
            }
        } // end construct local node element lists

        //- Communicate local elements to other processes
        vector<int> nCellsPerProc, cellNodeListSizes, nNodesPerProc;
        for(int proc = 0; proc < comm.nProcs(); ++proc)
        {
            nCellsPerProc.push_back(cellLists[proc].size());
            cellNodeListSizes.push_back(localCells[proc].size());
            nNodesPerProc.push_back(localNodes[proc].size());
        }

        printf("Sending local size info from main proc...\n");
        comm.scatter(comm.mainProcNo(), nCellsPerProc);
        comm.scatter(comm.mainProcNo(), cellNodeListSizes);
        comm.scatter(comm.mainProcNo(), nNodesPerProc);

        printf("Sending local mesh info from main proc...\n");
        for(int proc = 0; proc < comm.nProcs(); ++proc)
        {
            if(proc != comm.rank())
            {
                comm.send(proc, cellLists[proc], 0);
                comm.send(proc, localCellInds[proc], 1);
                comm.send(proc, localCells[proc], 2);
                comm.send(proc, localNodes[proc], 3);
            }
        }

        //- Send sendOrder/buffer info
        vector<int> nNeighbouringProcsPerProc(comm.nProcs());
        for(int i = 0; i < nNeighbouringProcsPerProc.size(); ++i)
            nNeighbouringProcsPerProc[i] = neighbouringProcs[i].size();

        //- Send neighbouring proc lists
        comm.scatter(comm.rank(), nNeighbouringProcsPerProc);
        for(int proc = 0; proc < comm.nProcs(); ++proc)
        {
            if(proc != comm.rank())
                comm.send(proc, neighbouringProcs[proc]);
        }

        for(int proc = 0; proc < comm.nProcs(); ++proc)
        {
            if(proc != comm.rank())
            {
                vector<Size> sendBufferSizes(neighbouringProcs[proc].size());
                vector<Size> recvBufferSizes(neighbouringProcs[proc].size());

                for(int i = 0; i < neighbouringProcs[proc].size(); ++i)
                {
                    int nbProc = neighbouringProcs[proc][i];
                    sendBufferSizes[i] = sendOrdering[proc][nbProc].size();
                    recvBufferSizes[i] = recvOrdering[proc][nbProc].size();
                }

                comm.send(proc, sendBufferSizes);
                comm.send(proc, recvBufferSizes);
            }
        }

        for(int proc = 0; proc < comm.nProcs(); ++proc)
        {
            if(proc != comm.rank())
            {
                for(int i = 0; i < neighbouringProcs[proc].size(); ++i)
                {
                    comm.send(proc, sendOrdering[proc][neighbouringProcs[proc][i]]);
                    comm.send(proc, recvOrdering[proc][neighbouringProcs[proc][i]]);
                }
            }
            else
            {
                for(int i = 0; i < neighbouringProcs[proc].size(); ++i)
                {
                    addNeighbouringProc(neighbouringProcs[proc][i],
                                        std::vector<Label>(),
                                        sendOrdering[proc][neighbouringProcs[proc][i]],
                            recvOrdering[proc][neighbouringProcs[proc][i]]);
                }
            }
        }

        comm.waitAll();
        init(localNodes[comm.rank()], localCellInds[comm.rank()], localCells[comm.rank()]);
    }
    else
    {
        //- General cell/mesh info
        int nCellsThisProc = comm.scatter(comm.mainProcNo(), vector<int>(1));
        int elementListSize = comm.scatter(comm.mainProcNo(), vector<int>(1));
        int nNodesThisProc = comm.scatter(comm.mainProcNo(), vector<int>(1));

        vector<Label> globalCellIds(nCellsThisProc);
        vector<Label> cellInds(nCellsThisProc + 1);
        vector<Label> cells(elementListSize);
        vector<Point2D> nodes(nNodesThisProc);

        comm.irecv(comm.mainProcNo(), globalCellIds, 0);
        comm.irecv(comm.mainProcNo(), cellInds, 1);
        comm.irecv(comm.mainProcNo(), cells, 2);
        comm.irecv(comm.mainProcNo(), nodes, 3);

        //- Get the number of partition neighbours this proc has
        int nNeighbouringProcs = comm.scatter(comm.mainProcNo(), vector<int>(1));

        //- Get a list of the neighbouring procs
        vector<int> neighbouringProcs(nNeighbouringProcs);
        comm.recv(comm.mainProcNo(), neighbouringProcs);

        //- Get the buffer size for each neighbouring proc
        vector<Size> sendBufferSizes(nNeighbouringProcs), recvBufferSizes(nNeighbouringProcs);
        comm.recv(comm.mainProcNo(), sendBufferSizes);
        comm.recv(comm.mainProcNo(), recvBufferSizes);

        for(int i = 0; i < nNeighbouringProcs; ++i)
        {
            vector<Label> sendOrder(sendBufferSizes[i]);
            vector<Label> recvOrder(recvBufferSizes[i]);

            comm.recv(comm.mainProc(), sendOrder);
            comm.recv(comm.mainProc(), recvOrder);

            addNeighbouringProc(neighbouringProcs[i], vector<Label>(), sendOrder, recvOrder);
        }

        comm.waitAll();
        init(nodes, cellInds, cells);
    }

    auto localSizes = comm.allGather(cells_.size());
    vector<int> globalInds(1, 0);
    for(int upper: localSizes)
        globalInds.push_back(globalInds.back() + upper);

    if(comm.mainProc())
    {
        for(int ind: globalInds)
            printf("%d\n", ind);
    }
}

void FiniteVolumeGrid2D::addNeighbouringProc(int procNo,
                                             const std::vector<Label> &bufferCells,
                                             const std::vector<Label> &sendOrder,
                                             const std::vector<Label> &recvOrder)
{
    neighbouringProcs_.push_back(procNo);
    //- Create a cell zone
    procSendOrder_.push_back(sendOrder);
    procRecvOrder_.push_back(recvOrder);
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

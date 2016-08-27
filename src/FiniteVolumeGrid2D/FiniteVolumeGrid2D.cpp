#include "FiniteVolumeGrid2D.h"
#include "Exception.h"

FiniteVolumeGrid2D::FiniteVolumeGrid2D(Size nNodes, Size nCells, Size nFaces)
    :
      nodeSearch_(nodes_)
{
    nodes_.reserve(nNodes);
    cells_.reserve(nCells);
    faces_.reserve(nFaces);
    fluidCells_.rename("FluidCells");
}

//- size info
std::string FiniteVolumeGrid2D::gridInfo() const
{
    using namespace std;

    return "Finite volume grid info:\n"
            "------------------------\n"
            "Number of nodes: " + to_string(nodes_.size()) + "\n"
            "Number of quad cells: " + to_string(quadCells_.size()) + "\n"
            "Number of tri cells: " + to_string(triCells_.size()) + "\n"
            "Number of cells total: " + to_string(cells_.size()) + "\n"
            "Number of interior faces: " + to_string(interiorFaces_.size()) + "\n"
            "Number of boundary faces: " + to_string(boundaryFaces_.size()) + "\n"
            "Number of faces total: " + to_string(faces_.size()) + "\n"
            "Bounding box: " + bBox_.toString() + "\n";
}

//- Create grid entities
Label FiniteVolumeGrid2D::createCell(const std::vector<Label>& nodeIds)
{
    using namespace std;

    cells_.push_back(Cell(nodeIds, nodes_));

    Cell &newCell = cells_.back();
    newCell.id_ = cells_.size() - 1;

    for(Label id: nodeIds)
        nodes_[id].cells_.push_back(std::cref(newCell));

    for(Label i = 0, end = nodeIds.size(); i < end; ++i)
    {
        Label n1 = nodeIds[i], n2 = nodeIds[(i + 1)%end];

        auto key = n1 < n2 ? make_pair(n1, n2) : make_pair(n2, n1);
        auto it = faceDirectory_.find(key);

        if(it == faceDirectory_.end()) // face doesn't exist, so create it
        {
            faces_.push_back(Face(n1, n2, nodes_, Face::BOUNDARY));

            Face &face = faces_.back();
            face.id_ = faces_.size() - 1;

            faceDirectory_[key] = face.id();

            face.addCell(newCell);
        }
        else // face already exists, but is now an interior face
        {
            Face &face = faces_[it->second];

            face.type_ = Face::INTERIOR;

            face.addCell(newCell);
        }
    }

    switch(nodeIds.size())
    {
    case 3:
        triCells_.push_back(std::cref(newCell));
        break;
    case 4:
        quadCells_.push_back(std::cref(newCell));
        break;
    };

    return newCell.id();
}

Label FiniteVolumeGrid2D::addNode(const Point2D &point)
{
    Label id = nodes_.size();
    nodes_.push_back(Node(point, id));
    return id;
}

//- Cell related methods

UniqueCellGroup& FiniteVolumeGrid2D::moveCellsToFluidCellGroup(const std::vector<size_t>& ids) const
{
    for(size_t id: ids)
    {
        cells_[id].setActive();
        cells_[id].setFluidCell();
        fluidCells_.moveToGroup(cells_[id]);
    }

    constructActiveCellGroup();
    return fluidCells_;
}

UniqueCellGroup& FiniteVolumeGrid2D::moveAllCellsToFluidCellGroup() const
{
    for(const Cell &cell: cells_)
    {
        cell.setActive();
        cell.setFluidCell();
    }

    fluidCells_.moveAllCellsToThisGroup();

    constructActiveCellGroup();
    return fluidCells_;
}

UniqueCellGroup& FiniteVolumeGrid2D::moveCellsToInactiveCellGroup(const std::vector<size_t>& ids) const
{
    for(size_t id: ids)
    {
        cells_[id].setInactive();
        cells_[id].setNonFluidCell();
        inactiveCells_.moveToGroup(cells_[id]);
    }

    constructActiveCellGroup();
    return inactiveCells_;
}

UniqueCellGroup& FiniteVolumeGrid2D::moveCellsToCellGroup(const std::string& name, const std::vector<size_t>& ids) const
{
    UniqueCellGroup &group = (cellGroups_.insert(std::make_pair(name, UniqueCellGroup(name)))).first->second;

    for(size_t id: ids)
    {
        cells_[id].setActive();
        cells_[id].setNonFluidCell();
        group.moveToGroup(cells_[id]);
    }

    constructActiveCellGroup();
    return group;
}

CellGroup& FiniteVolumeGrid2D::cellGroup(const std::string &name)
{
    if(name == "fluidCells")
        return fluidCells_;
    else if(name == "activeCells")
        return activeCells_;
    else
        return cellGroups_[name];
}

const CellGroup& FiniteVolumeGrid2D::cellGroup(const std::string &name) const
{
    if(name == "fluidCells")
        return fluidCells_;
    else if(name == "activeCells")
        return activeCells_;
    else
        return cellGroups_.find(name)->second;
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

void FiniteVolumeGrid2D::applyPatch(const std::string &patchName, const std::vector<Ref<Face> > &faces)
{
    auto insert = patches_.insert(std::make_pair(patchName, Patch(patchName, patches_.size())));

    if(!insert.second)
        throw Exception("FiniteVolumeGrid2D", "applyPatch", "patch already exists.");

    Patch &patch = (insert.first)->second;

    for(Face &face: faces)
        patch.addFace(face);
}

const Node& FiniteVolumeGrid2D::findNearestNode(const Point2D& pt) const
{
    return nodeSearch_.kNearestNeighbourSearch(pt, 1)[0];
}

//- Protected methods

void FiniteVolumeGrid2D::initNodes()
{
    nodeSearch_.constructRTree();
}

void FiniteVolumeGrid2D::initCells()
{  
    for(Cell &cell: cells_)
    {
        cell.isActive_ = true;
        fluidCells_.push_back(cell);
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

    constructActiveCellGroup();
}

void FiniteVolumeGrid2D::initConnectivity()
{
    initNodes();
    initCells();
}

void FiniteVolumeGrid2D::constructActiveCellGroup() const
{
    Index idx = 0;
    activeCells_.clear();
    activeCells_.reserve(cells_.size());

    for(const Cell &cell: cells_)
    {
        if(cell.isActive())
        {
            activeCells_.push_back(cell);
            cell.globalIndex_ = idx++; // compute the global index
        }
        else
            cell.globalIndex_ = Cell::INACTIVE;
    }
}

void FiniteVolumeGrid2D::computeBoundingBox()
{
    bBox_ = BoundingBox(nodes_.data(), nodes_.size());
}

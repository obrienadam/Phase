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
            "Number of cells: " + to_string(cells_.size()) + "\n"
            "Number of interior faces: " + to_string(interiorFaces_.size()) + "\n"
            "Number of boundary faces: " + to_string(boundaryFaces_.size()) + "\n"
            "Number of faces total: " + to_string(faces_.size()) + "\n"
            "Bounding box: " + bBox_.toString() + "\n";
}

//- Create grid entities
Label FiniteVolumeGrid2D::createCell(const std::vector<Label>& nodeIds)
{
    using namespace std;

    Label id = cells_.size();

    cells_.push_back(Cell(nodeIds, nodes_));

    Cell &newCell = cells_.back();
    newCell.id_ = id;

    for(Label i = 0, end = nodeIds.size(); i < end; ++i)
    {
        Label n1 = nodeIds[i], n2 = nodeIds[(i + 1)%end];

        auto key = n1 < n2 ? make_pair(n1, n2) : make_pair(n2, n1);
        auto it = faceDirectory_.find(key);

        if(it == faceDirectory_.end()) // face doesn't exist, so create it
        {
            Face face = Face(n1, n2, nodes_, Face::BOUNDARY);
            face.id_ = faces_.size();

            faces_.push_back(face);
            faceDirectory_[key] = face.id();

            faces_.back().addCell(newCell);
            newCell.boundaryLinks_.push_back(BoundaryLink(newCell, faces_.back()));
        }
        else // face already exists, but is now an interior face
        {
            Face &face = faces_[it->second];

            face.type_ = Face::INTERIOR;

            Cell &cell = cells_[face.lCell().id()];

            for(auto it = cell.boundaryLinks_.begin(); it != cell.boundaryLinks_.end(); ++it)
            {
                if(it->face().id() == face.id())
                {
                    cell.boundaryLinks_.erase(it);
                    break;
                }
            }

            face.addCell(newCell);

            cell.interiorLinks_.push_back(InteriorLink(cell, face, newCell));
            newCell.interiorLinks_.push_back(InteriorLink(newCell, face, cell));
        }
    }

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
        cells_[id].isActive_ = true;
        fluidCells_.moveToGroup(cells_[id]);
    }

    constructActiveCellGroup();
    return fluidCells_;
}

UniqueCellGroup& FiniteVolumeGrid2D::moveAllCellsToFluidCellGroup() const
{
    for(const Cell &cell: cells_)
        cell.isActive_ = true;

    fluidCells_.moveAllCellsToThisGroup();

    constructActiveCellGroup();
    return fluidCells_;
}

UniqueCellGroup& FiniteVolumeGrid2D::moveCellsToInactiveCellGroup(const std::vector<size_t>& ids) const
{
    for(size_t id: ids)
    {
        cells_[id].isActive_ = false;
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
        cells_[id].isActive_ = true;
        group.moveToGroup(cells_[id]);
    }

    constructActiveCellGroup();
    return group;
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

    return it->second;
}

void FiniteVolumeGrid2D::applyPatch(const std::string &patchName, const std::vector<Ref<Face> > &faces)
{
    auto insert = patches_.insert(std::make_pair(patchName, Patch(patches_.size(), patchName)));

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
            boundaryFaces_.push_back(face);
        else
            interiorFaces_.push_back(face);
    }

    constructActiveCellGroup();
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

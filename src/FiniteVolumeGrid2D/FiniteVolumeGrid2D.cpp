#include "FiniteVolumeGrid2D.h"
#include "Exception.h"

FiniteVolumeGrid2D::FiniteVolumeGrid2D(size_t nNodes, size_t nCells, size_t nFaces)
{
    nodes_.reserve(nNodes);
    cells_.reserve(nCells);
    faces_.reserve(nFaces);
    patches_.reserve(4);
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
size_t FiniteVolumeGrid2D::createFace(size_t lNodeId, size_t rNodeId, Face::Type type)
{
    using namespace std;
    typedef pair<size_t, size_t> Key;

    Face face(lNodeId, rNodeId, nodes_, type);
    face.id_ = faces_.size();

    faces_.push_back(face);

    Key key1(lNodeId, rNodeId), key2(rNodeId, lNodeId);

    faceDirectory_.insert(pair<Key, size_t>(key1, face.id()));
    faceDirectory_.insert(pair<Key, size_t>(key2, face.id()));

    if(type == Face::INTERIOR)
        interiorFaces_.push_back(Ref<const Face>(faces_.back()));
    else if (type == Face::BOUNDARY)
        boundaryFaces_.push_back(Ref<const Face>(faces_.back()));

    return face.id();
}

size_t FiniteVolumeGrid2D::createCell(const std::vector<size_t> &faceIds)
{
    Cell cell(faceIds, faces_);
    cell.id_ = cells_.size();
    cell.globalIndex_ = cell.id();
    cells_.push_back(cell);

    for(size_t faceId: faceIds)
        faces_[faceId].addCell(cells_.back());

    return cell.id();
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

void FiniteVolumeGrid2D::applyPatch(const std::string &patchName, const std::vector<Ref<Face> > &faces)
{
    patches_.push_back(Patch(patches_.size(), patchName));

    for(Face &face: faces)
        patches_.back().addFace(face);
}

void FiniteVolumeGrid2D::computeBoundingBox()
{
    bBox_ = BoundingBox(nodes_.data(), nodes_.size());
}

//- Protected methods

void FiniteVolumeGrid2D::initCells()
{
    for(Cell &cell: cells_)
    {
        cell.computeCellAdjacency();
        cell.isActive_ = true;
        fluidCells_.push_back(cell);
    }

    constructActiveCellGroup();
}

void FiniteVolumeGrid2D::constructActiveCellGroup() const
{
    size_t idx = 0;
    activeCells_.clear();
    activeCells_.reserve(cells_.size());

    for(const Cell &cell: cells_)
    {
        if(cell.isActive())
        {
            activeCells_.push_back(cell);
            cell.globalIndex_ = idx++; // compute the global index
        }
    }
}

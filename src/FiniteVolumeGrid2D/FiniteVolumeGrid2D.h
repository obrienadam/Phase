#ifndef FINITE_VOLUME_GRID_2D_H
#define FINITE_VOLUME_GRID_2D_H

#include <vector>
#include <utility>
#include <map>

#include "Node.h"
#include "Cell.h"
#include "Face.h"
#include "Patch.h"
#include "BoundingBox.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

class FiniteVolumeGrid2D
{
public:

    FiniteVolumeGrid2D(size_t nNodes = 0, size_t nCells = 0, size_t nFaces = 0);

    size_t nNodes() const { return nodes.size(); }
    size_t nCells() const { return cells.size(); }
    size_t nFaces() const { return faces.size(); }

    size_t nActiveCells() const { return cells.size(); }

    size_t createFace(size_t lNodeId, size_t rNodeId, Face::Type type = Face::INTERIOR);
    size_t createCell(const std::vector<size_t>& faceIds);

    bool faceExists(size_t lNodeId, size_t rNodeId) const;
    size_t findFace(size_t lNodeId, size_t rNodeId) const;

    std::string gridInfo() const;

    void computeCellAdjacency();
    size_t computeGlobalIndices();

    void addPatch(const std::string& patchName);
    const std::vector<Patch>& patches() const { return patches_; }

    ScalarFiniteVolumeField& addScalarField(const std::string& fieldName) const;
    VectorFiniteVolumeField& addVectorField(const std::string& fieldName) const;

    ScalarFiniteVolumeField& addScalarField(const Input& input, const std::string& fieldName) const;
    VectorFiniteVolumeField& addVectorField(const Input& input, const std::string& fieldName) const;

    std::map<std::string, ScalarFiniteVolumeField>& scalarFields() const { return scalarFields_; }
    std::map<std::string, VectorFiniteVolumeField>& vectorFields() const { return vectorFields_; }

    std::vector<Node> nodes;
    std::vector<Cell> cells;
    std::vector<Face> faces;

protected:

    void computeBoundingBox();
    void applyPatch(const std::string& patchName, const std::vector< Ref<Face> >& faces);

    std::map<std::pair<size_t, size_t>, size_t> faceDirectory_;

    std::vector<size_t> interiorFaces_;
    std::vector<size_t> boundaryFaces_;

    std::vector<Patch> patches_;

    mutable std::map<std::string, ScalarFiniteVolumeField> scalarFields_;
    mutable std::map<std::string, VectorFiniteVolumeField> vectorFields_;

    BoundingBox bBox_;
};

#endif

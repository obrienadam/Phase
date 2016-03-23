#include "RectilinearGrid2D.h"
#include "Patch.h"

RectilinearGrid2D::RectilinearGrid2D(int nCellsI, int nCellsJ, Scalar hx, Scalar hy)
    :
      FiniteVolumeGrid2D((nCellsI + 1)*(nCellsJ + 1), nCellsI*nCellsJ, (nCellsI + 1)*(nCellsJ) + (nCellsI)*(nCellsJ + 1))
{
    hx_ = hx;
    hy_ = hy;

    int nNodesI = nCellsI + 1, nNodesJ = nCellsJ + 1;
    int uNodeI = nNodesI - 1, uNodeJ = nNodesJ - 1;

    auto nodeId = [nNodesI](int i, int j){ return i + j*nNodesI; };

    for(int j = 0; j < nNodesJ; ++j)
        for(int i = 0; i < nNodesI; ++i)
        {
            nodes.push_back(Node(i*hx_, j*hy_, nodeId(i, j)));

            if(i > 0)
            {
                createFace(nodeId(i - 1, j), nodeId(i, j), j == 0 || j == uNodeJ ? Face::BOUNDARY : Face::INTERIOR);

                if(j == 0)
                    bottomFaces_.push_back(faces.back());

                else if(j == uNodeJ)
                    topFaces_.push_back(faces.back());
            }

            if(j > 0)
            {
                createFace(nodeId(i, j - 1), nodeId(i, j), i == 0 || i == uNodeI ? Face::BOUNDARY : Face::INTERIOR);

                if(i == 0)
                    leftFaces_.push_back(faces.back());

                else if(i == uNodeI)
                    rightFaces_.push_back(faces.back());
            }
        }

    std::vector<size_t> faceIds(4);
    for(int j = 0; j < nCellsJ; ++j)
        for(int i = 0; i < nCellsI; ++i)
        {   
            faceIds[0] = findFace(nodeId(i + 1, j), nodeId(i + 1, j + 1));
            faceIds[1] = findFace(nodeId(i + 1, j + 1), nodeId(i, j + 1));
            faceIds[2] = findFace(nodeId(i, j + 1), nodeId(i, j));
            faceIds[3] = findFace(nodeId(i, j), nodeId(i + 1, j));

            createCell(faceIds);
        }

    applyBottomPatch("y-");
    applyTopPatch("y+");
    applyLeftPatch("x-");
    applyRightPatch("x+");

    computeCellAdjacency();
    computeBoundingBox();
}

void RectilinearGrid2D::applyBottomPatch(const std::string &patchName)
{
    applyPatch(patchName, bottomFaces_);
}

void RectilinearGrid2D::applyTopPatch(const std::string &patchName)
{
    applyPatch(patchName, topFaces_);
}

void RectilinearGrid2D::applyLeftPatch(const std::string& patchName)
{
    applyPatch(patchName, leftFaces_);
}

void RectilinearGrid2D::applyRightPatch(const std::string &patchName)
{
    applyPatch(patchName, rightFaces_);
}



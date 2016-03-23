#ifndef RECTILINEAR_GRID_2D
#define RECTILINEAR_GRID_2D

#include "FiniteVolumeGrid2D.h"

class RectilinearGrid2D : public FiniteVolumeGrid2D
{
public:

    RectilinearGrid2D(int nCellsI, int nCellsJ, Scalar hx, Scalar hy);

    void applyBottomPatch(const std::string& patchName);
    void applyTopPatch(const std::string& patchName);
    void applyLeftPatch(const std::string& patchName);
    void applyRightPatch(const std::string& patchName);

private:



    Scalar hx_, hy_;
    std::vector< Ref<Face> > bottomFaces_, topFaces_, leftFaces_, rightFaces_;
};

#endif

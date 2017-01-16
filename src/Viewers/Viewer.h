#ifndef VIEWER_H
#define VIEWER_H

#include <vector>

#include "Solver.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

class Viewer
{
public:

    Viewer(const Solver& solver, const Input& input, const std::string &customName = "");

    virtual void write(Scalar solutionTime);
    virtual void write(const std::vector<VolumeIntegrator>& volumeIntegrators);
    void addPolygons(const std::vector<Polygon>& polygons, const std::string& zoneName = "Polygon Zone");

protected:

    int  createZone(const FiniteVolumeGrid2D& grid);
    void writeCoords(const FiniteVolumeGrid2D& grid);
    int  writeConnectivity(const FiniteVolumeGrid2D& grid);
    void writeBoundaryConnectivity(const FiniteVolumeGrid2D& grid);
    void writeImmersedBoundaries(const Solver& solver);
    int  updateFlowSolutionPointers(Scalar solutionTime);

    const Solver& solver_;

    int fileId_, baseId_, zoneId_;
    std::vector<Scalar> timeValues_;
    std::vector<std::string> flowSolutionPointers_;

    std::vector< Ref<const ScalarFiniteVolumeField> > scalarFields_;
    std::vector< Ref<const VectorFiniteVolumeField> > vectorFields_;

    std::string caseName_, outputFilename_;
};

#endif
